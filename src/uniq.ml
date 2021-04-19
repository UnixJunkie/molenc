
(* uniq filter: keep only line if given field was never seen before  *)

open Printf

module CLI = Minicli.CLI
module Db = Dokeysto_camltc.Db_camltc.RW
module Ht = Hashtbl
module Log = Dolog.Log
module String = BatString
module Utls = Molenc.Utls

module type HT = sig
  type t
  val create: string -> t
  val mem: t -> string -> bool
  val add: t -> string -> unit
  val close: t -> unit
end

module HtOnDisk: HT = struct

  type t = Dokeysto_camltc.Db_camltc.RW.t
  
  let create input_fn =
    Db.create (input_fn ^ ".uniq.db")

  let mem db field =
    Db.mem db field

  let add db field =
    Db.add db field ""

  let close db =
    Db.close db
end

module HtInRAM: HT = struct

  type t = (string, unit) Ht.t
  
  let create input_fn =
    Ht.create (Utls.count_lines_of_file input_fn)

  let mem db field =
    Ht.mem db field

  let add db field =
    Ht.add db field ()

  let close _db =
    ()
end

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
               -f <int>: field to filter on\n  \
               [--in-RAM]: Ht in RAM rather than on disk\n"
        Sys.argv.(0);
      exit 1
    end;
  let mod_db =
    if CLI.get_set_bool ["--in-RAM"] args then
      (module HtInRAM: HT)
    else (module HtOnDisk: HT) in
  let module DB = (val mod_db: HT) in
  let input_fn = CLI.get_string ["-i"] args in
  let db_fn = input_fn ^ ".uniq.db" in
  let db = DB.create db_fn in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field = (CLI.get_int ["-f"] args) - 1 in
  Utls.with_in_file input_fn (fun input ->
      try
        let count = ref 0 in
        while true do
          let line = input_line input in
          let field = String.cut_on_char sep field line in
          (if not (DB.mem db field) then
             (DB.add db field;
              printf "%s\n" line)
          );
          incr count;
          (if !count mod 1000 = 0 then
             eprintf "done: %d\r%!" !count);
        done
      with End_of_file -> DB.close db
    )

let () = main ()
