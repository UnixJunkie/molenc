(* Efficiently extract tags from .sdf[.gz] files (sdf2csv) *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString
module StringSet = BatSet.String

let sdf_tag_regexp = Str.regexp "^>  <\\(.+\\)>"

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.sdf[.gz] -o molecules.csv -t tag1[,tag2[,...]\n  \
              -i <filename>: SDF (or .sdf.gz) input file\n  \
              -o <filename>: CSV output file\n  \
              -t <string>: comma-separated list of tags to extract\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  Log.set_log_level (if verbose then Log.DEBUG else Log.INFO);
  Log.color_on ();
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let tags = CLI.get_string ["-t"] args in
  (* ----------------------------------------------------------------------- *)
  let tags_list = S.split_on_char ',' tags in
  let tags_set = StringSet.of_list tags_list in
  let num_tags = L.length tags_list in
  (* FBR:TODO compressed input *)
  let tag2values = Ht.create num_tags in
  L.iter (fun tag ->
      Ht.add tag2values tag []
    ) tags_list;
  let input = open_in input_fn in
  (try
     while true do
       let line = input_line input in
       if S.length line > 0 && String.unsafe_get line 0 = '>' then
         let _i = Str.search_forward sdf_tag_regexp line 0 in
         let tag = Str.matched_group 1 line in
         if StringSet.mem tag tags_set then
           let prev = Ht.find tag2values tag in
           (* the next line in the SDF has the actual value for that tag *)
           let value = input_line input in
           Ht.add tag2values tag (value :: prev)
     done
   with End_of_file -> close_in input
  );
  (* remove pipe file? *)
  let expected = match tags_list with
    | [] -> failwith "Get_sdf_tag.main: empty list of tags"
    | tag0 :: _tags -> L.length (Ht.find tag2values tag0) in
  let values = A.make num_tags [|""|] in
  L.iteri (fun i tag ->
      let vals = Ht.find tag2values tag in
      let arr = A.of_list vals in
      let n = A.length arr in
      (* check all tags were seen the same number of times *)
      (if n <> expected then
         (Log.fatal "expected: %d <> len(res[%s]) = %d" expected tag n;
          exit 1)
      );
      (* values need to be reversed because they were accumulated in a list *)
      A.rev_in_place arr;
      values.(i) <- arr
    ) tags_list;
  (* output CSV file *)
  LO.with_out_file output_fn (fun output ->
      for i = 0 to expected - 1 do
        for j = 0 to num_tags - 1 do
          fprintf output "%s," values.(j).(i)
        done;
        output_char output '\n' (* EOL *)
      done
    )

let () = main ()
