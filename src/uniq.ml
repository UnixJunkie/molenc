
(* uniq filter: keep only line if given field was never seen before  *)

open Printf

module CLI = Minicli.CLI
module Ht = Hashtbl
module Log = Dolog.Log
module String = BatString
module Utls = Molenc.Utls

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
               -f <int>: field to filter on\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field = (CLI.get_int ["-f"] args) - 1 in
  Utls.with_in_file input_fn (fun input ->
      try
        let ht = Hashtbl.create 11 in
        while true do
          let line = input_line input in
          let field = String.cut_on_char sep field line in
          if not (Ht.mem ht field) then
            begin
              Ht.add ht field ();
              printf "%s\n" line
            end
        done
      with End_of_file -> ()
    )

let () = main ()
