
(* Compute the Pareto front for a set of solutions *)

open Printf

module CLI = Minicli.CLI
module L = BatList
module Log = Dolog.Log

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
               -f <int>,<int>[,...]: fields to filter on\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field_num = (CLI.get_int ["-f"] args) - 1 in
  (* FBR: skip lines starting with a '#' *)
  failwith "not implemented yet"

let () = main ()
