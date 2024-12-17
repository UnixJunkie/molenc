(* replacement for UNIX's shuf command, but seedable for reproducibility *)

open Printf

module CLI = Minicli.CLI
module Log = Dolog.Log
module LO = Line_oriented

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file\n  \
               -o <filename>: output file\n  \
               [-n <int>]: output at most N lines (default=all)\n  \
               [-s <int>]: random seed (default=none)\n"
        Sys.argv.(0);
      exit 1
    end;
  CLI.finalize ();
  let sorted = CLI.get_set_bool ["--sorted"] args in
  let force = CLI.get_set_bool ["--force"] args in
  let input_fn = CLI.get_string ["-i"] args in
  failwith "not implemented yet"

let () = main ()
