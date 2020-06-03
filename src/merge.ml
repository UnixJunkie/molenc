
(* merge score files using average of normalized scores
   (unweighted score merging strategy) *)

open Printf

module CLI = Minicli.CLI
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module String = BatString
module Utls = Molenc.Utls

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -ifs <filename:sep:name_field:{-}score_field>,...: input score files\n  \
              the optional '-' means that lower scores are better (like docking scores)\n  \
              -o <filename>: output scores file\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize();
  let all_scores = ref [] in
  (* read all scores *)
  Utls.iter_on_lines_of_file input_fn (fun line ->
      let score_field = String.cut_on_char sep field line in
      let score =
        try Scanf.sscanf score_field "%f" (fun x -> x)
        with exn ->
          begin
            Log.fatal "Rank: cannot parse float: %s" score_field;
            raise exn
          end in
      all_scores := score :: !all_scores
    );
  failwith "not implemented yet"

let () = main ()
