
(* uniq filter: keep only line if given field was never seen before  *)

open Printf

module CLI = Minicli.CLI
module Db = Dokeysto.Db.RW
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
  let db_fn = input_fn ^ ".uniq.db" in
  let db = Db.create db_fn in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field = (CLI.get_int ["-f"] args) - 1 in
  Utls.with_in_file input_fn (fun input ->
      try
        let count = ref 0 in
        while true do
          let line = input_line input in
          let field = String.cut_on_char sep field line in
          (if not (Db.mem db field) then
             (Db.add db field "";
              printf "%s\n" line)
          );
          incr count;
          (if !count mod 1000 = 0 then
             eprintf "done: %d\r%!" !count);
        done
      with End_of_file -> Db.close db
    )

let () = main ()
