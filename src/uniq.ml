
(* uniq filter: keep only line if given field was never seen before  *)

open Printf

module CLI = Minicli.CLI
module Db = Dokeysto.Db.RW
module Ht = Hashtbl
module Log = Dolog.Log
module LO = Line_oriented
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

  type t = Dokeysto.Db.RW.t

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

  let create _input_fn =
    (* DO NOT ever try to read the whole input file in case of --in-RAM:
     * we want to be able to read molecules from a UNIX pipe *)
    Ht.create 1_000_000

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
               [--force]: erase index files, if any\n  \
               [--sorted]: file already sorted on that field\n  \
               [--in-RAM]: Ht in RAM rather than on disk\n"
        Sys.argv.(0);
      exit 1
    end;
  let mod_db =
    if CLI.get_set_bool ["--in-RAM"] args then
      (module HtInRAM: HT)
    else (module HtOnDisk: HT) in
  let module DB = (val mod_db: HT) in
  let sorted = CLI.get_set_bool ["--sorted"] args in
  let force = CLI.get_set_bool ["--force"] args in
  let input_fn = CLI.get_string ["-i"] args in
  (if force then
     (Utls.rm_file (input_fn ^ ".uniq.db");
      Utls.rm_file (input_fn ^ ".uniq.db.idx"))
  );
  let db = DB.create input_fn in
  let prev_field = ref "" in
  let uniq_field_check, register_field =
    if sorted then
      ((fun field -> !prev_field <> field),
       (fun field -> prev_field := field))
    else
      ((fun field -> not (DB.mem db field)),
       (fun field -> DB.add db field))
  in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field_num = (CLI.get_int ["-f"] args) - 1 in
  let count = ref 0 in
  LO.with_in_file input_fn (fun input ->
      try
        while true do
          let line = input_line input in
          let field_str = String.cut_on_char sep field_num line in
          (if uniq_field_check field_str then
             (register_field field_str;
              printf "%s\n" line)
          );
          incr count;
          (if !count mod 1000 = 0 then
             eprintf "done: %d\r%!" !count
          )
        done
      with End_of_file -> DB.close db
    )

let () = main ()
