
(* Compute the Pareto front for a set of solutions *)

open Printf

module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

let parse_field_nums s =
  let strings = S.split_on_char ',' s in
  (* the user wants 1st field at index 1 (as in Awk), but on computers
     the 1st index is 0 *)
  L.map (fun num_str -> (int_of_string num_str) - 1) strings

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file\n  \
               -d <char>: field separator (default=\\t)\n  \
               -f <int>,<int>[,...]: fields to filter on \
               (first field index=1; same as awk)\n" Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  (* ignore comment lines (starting with a '#') *)
  let _all_lines = LO.filter input_fn (fun l -> not (S.starts_with l "#")) in
  let _sep = CLI.get_char_def ["-d"] args '\t' in
  let field_nums_str = CLI.get_string ["-f"] args in
  CLI.finalize();
  let _field_nums = parse_field_nums field_nums_str in
  failwith "not implemented yet"

let () = main ()
