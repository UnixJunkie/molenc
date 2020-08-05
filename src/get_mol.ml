(* extract molecules with given names from a MOL2, SDF or SMILES file
   molecules order is preserved and follows the one of provided names *)

open Printf

module CLI = Minicli.CLI
module DB = Dokeysto_camltc.Db_camltc.RW
module L = BatList
module Log = Dolog.Log
module S = BatString
module Utls = Molenc.Utls

let db_name_of fn =
  fn ^ ".db"

type mol_name_provider = On_cli of string
                       | From_file of string

let mol_reader_for_file fn =
  if S.ends_with fn ".mol2" then Mol2.(read_one_raw, get_name)
  else if S.ends_with fn ".sdf" then Sdf.(read_one, get_fst_line)
  else if S.ends_with fn ".smi" then Smi.(read_one, get_name)
  else failwith ("Get_mol.mol_reader_for_file: not {.mol2|.sdf|.smi}: " ^ fn)

let populate_db db input_fn =
  let read_one_mol, read_mol_name = mol_reader_for_file input_fn in
  let count = ref 0 in
  Utls.with_in_file input_fn (fun input ->
      try
        while true do
          let m = read_one_mol input in
          Log.debug "m: %s" m;
          let name = read_mol_name m in
          Log.debug "name: %s" name;
          DB.add db name m;
          incr count;
          if (!count mod 10_000) = 0 then
            eprintf "read %d\r%!" !count;
        done
      with End_of_file -> DB.sync db
    )

let db_open_or_create verbose force input_fn =
  let db_fn = db_name_of input_fn in
  (* is there a DB already? *)
  let db_exists, db =
    if force || not (Sys.file_exists db_fn) then
      (Log.info "creating %s" db_fn;
       (false, DB.create db_fn))
    else
      (Log.warn "reusing %s" db_fn;
       (true, DB.open_existing db_fn)) in
  if not db_exists then populate_db db input_fn;
  if verbose then
    DB.iter (fun k v ->
        Log.debug "k: %s v: %s" k v
      ) db;
  db

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.{sdf|mol2|smi} \
              {-names \"mol1,mol2,...\"|-f names_file} [-v]\n  \
              -i <filename>: molecules input file\n  \
              [-o <filename>]: molecules output file (default=stdout)\n  \
              [-names <string>,<string>,...]: molecule names\n  \
              [-f <filename>]: get molecule names from file\n  \
              [-if <filename>,<filename>,...]: several molecule input files\n  \
              [--force]: overwrite existing db file(s), if any\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  Log.set_log_level (if verbose then Log.DEBUG else Log.INFO);
  Log.set_output stderr;
  Log.color_on ();
  let input_fns =
    match (CLI.get_string_opt ["-i"] args, CLI.get_string_opt ["-if"] args) with
    | (None, None) -> failwith "Get_mol: provide either -i or -if"
    | (Some _, Some _) -> failwith "Get_mol: both -i and -if"
    | (Some fn, None) -> [fn]
    | (None, Some filenames) -> S.split_on_string filenames ~by:"," in
  let maybe_output_fn = CLI.get_string_opt ["-o"] args in
  let force_db_creation = CLI.get_set_bool ["--force"] args in
  let names_provider = match CLI.get_string_opt ["-names"] args with
    | Some names -> On_cli names
    | None ->
      let fn = CLI.get_string ["-f"] args in
      From_file fn in
  CLI.finalize ();
  let names = match names_provider with
    | On_cli names -> S.split_on_string names ~by:","
    | From_file fn -> Utls.lines_of_file fn in
  let dbs = L.map (db_open_or_create verbose force_db_creation) input_fns in
  let out = match maybe_output_fn with
    | None -> stdout
    | Some output_fn -> open_out_bin output_fn in
  List.iter (fun name ->
      try
        (* find containing db, if any *)
        let db = L.find (fun db -> DB.mem db name) dbs in
        (* extract molecule from it *)
        let m = DB.find db name in
        fprintf out "%s" m
      with Not_found ->
        (* no db contains this molecule *)
        Log.warn "not found: %s" name
    ) names;
  L.iter DB.close dbs;
  (match maybe_output_fn with
   | Some _fn -> close_out out
   | None -> ())

let () = main ()
