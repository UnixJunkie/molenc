(* extract molecules with given names from a MOL2, SDF or SMILES file
   molecules order is preserved and follows the one of provided names *)

open Printf

module CLI = Minicli.CLI
module DB = Dokeysto.Db.RW
module Fn = Filename
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString
module StringSet = BatSet.String
module Utls = Molenc.Utls

type mol_name_provider = On_cli of string
                       | From_file of string

let mol_reader_for_file fn =
  if S.ends_with fn ".mol2" || S.ends_with fn ".mol2.gz" then
    Mol2.(read_one_raw, get_name)
  else if S.ends_with fn ".sdf" || S.ends_with fn ".sdf.gz" then
    Sdf.(read_one, get_fst_line)
  else if S.ends_with fn ".smi" || S.ends_with fn ".smi.gz" then
    Smi.(read_one, get_name)
  else if S.ends_with fn ".ph4" || S.ends_with fn ".ph4.gz" then
    Ph4.(read_one, get_name)
  else failwith ("Get_mol.mol_reader_for_file: not \
                  {.mol2[.gz]|.sdf[.gz]|.smi[.gz]|.ph4[.gz]}: " ^ fn)

(* like LO.with_in_file but transparently supporting .gz files *)
let with_in_file fn f =
  let input =
    if S.ends_with fn ".gz" then
      let cmd = sprintf "zcat %s" fn in
      Log.info "running: %s" cmd;
      Unix.open_process_in cmd
    else
      open_in_bin fn in
  let res = f input in
  close_in input;
  res

(* several consecutive molecules w/ same name as a large string *)
let cat_mols mols =
  String.concat "" (L.rev mols)

(* read input channel in full *)
let string_of_in_chan c =
  really_input_string c (in_channel_length c)

let lz4_string tmp_in_fn tmp_out_fn compress_or_not s =
  (* dump string to file *)
  LO.with_out_file tmp_in_fn (fun out ->
      output_string out s
    );
  (* compress/decompress *)
  let cmd =
    sprintf (if compress_or_not
             then "lz4 -zcq %s > %s"
             else "lz4 -dcq %s > %s")
      tmp_in_fn tmp_out_fn in
  (* Log.info "running: %s" cmd; *)
  let ret = Sys.command cmd in
  if ret <> 0 then
    (Log.fatal "Get_mol.lz4_string: command %s terminated w/ %d"
       cmd ret;
     exit 1)
  else
    (* read lz4 output *)
    LO.with_in_file tmp_out_fn string_of_in_chan

(* consecutive molecules w/ same name are stored under the same entry *)
let read_all_molecules compress db_add db_close input_fn =
  let read_one_mol, read_mol_name = mol_reader_for_file input_fn in
  let count = ref 0 in
  let mols = ref [] in
  let prev_name = ref "" in
  let tmp_in_fn, tmp_out_fn, lz4 =
    if compress then
      let tmp_in_fn = Fn.temp_file ~temp_dir:"/tmp" "get_mol_" "" in
      let tmp_out_fn = Fn.temp_file ~temp_dir:"/tmp" "get_mol_" ".lz4" in
      (tmp_in_fn, tmp_out_fn, lz4_string tmp_in_fn tmp_out_fn true)
    else
      ("", "", fun x -> x) in
  with_in_file input_fn (fun input ->
      try
        while true do
          let m = read_one_mol input in
          Log.debug "m: %s" m;
          let name = read_mol_name m in
          Log.debug "name: %s" name;
          (if name <> !prev_name then
             begin
               if !mols <> [] then
                 db_add !prev_name (lz4 (cat_mols !mols))
               ;
               mols := [m];
               prev_name := name
             end
           else
             mols := m :: !mols
          );
          incr count;
          if (!count mod 10_000) = 0 then
            eprintf "read %d\r%!" !count;
        done;
        if !mols <> [] then
          db_add !prev_name (lz4 (cat_mols !mols))
      with End_of_file -> db_close ()
    );
  if tmp_in_fn <> "" then Unix.unlink tmp_in_fn;
  if tmp_out_fn <> "" then Unix.unlink tmp_out_fn

let populate_db compress db input_fn =
  let db_add name m =
    DB.add db name m in
  let db_close () =
    DB.sync db in
  read_all_molecules compress db_add db_close input_fn

(* almost copy/paste of populate_db above ... *)
let populate_ht names input_fn =
  let required_names = StringSet.of_list names in
  let nb_names = StringSet.cardinal required_names in
  let collected = Ht.create nb_names in
  let db_add name m =
    if StringSet.mem name required_names then
      Ht.add collected name m in
  let db_close () = () in
  (* for an in-RAM ht, never compress *)
  read_all_molecules false db_add db_close input_fn;
  collected

let db_open_or_create verbose force compress input_fn =
  let db_fn =
    if compress then
      input_fn ^ ".lz4db"
    else
      input_fn ^ ".db" in
  (* is there a DB already? *)
  let db_exists, db =
    if force || not (Sys.file_exists db_fn) then
      (Log.info "creating %s" db_fn;
       Utls.rm_file db_fn;
       Utls.rm_file (db_fn ^ ".idx");
       (false, DB.create db_fn))
    else
      (Log.warn "reusing %s" db_fn;
       (true, DB.open_existing db_fn)) in
  if not db_exists then populate_db compress db input_fn;
  if verbose then
    DB.iter (fun k v ->
        Log.debug "k: %s v: %s" k v
      ) db;
  db

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.{sdf|mol2|smi|ph4}[.gz] \
              {-names \"mol1,mol2,...\"|-f names_file} [-v]\n  \
              -i <filename>: molecules input file\n  \
              [-o <filename[.gz]>]: molecules output file (default=stdout)\n  \
              [-names <string>,<string>,...]: molecule names\n  \
              [-f <filename>]: get molecule names from file\n  \
              [-nf <filename>]: output unfound molecule names to file\n  \
              [-if <filename>]: read molecule input file names from file\n  \
              [--force]: overwrite existing db file(s), if any\n  \
              [--no-index]: do not create index files to accelerate\n  \
              future queries on same input files\n\  \
              [-z]: toggle unlz4/lz4 compression of values from/to cache files\n"
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
    | (None, Some fn) -> LO.lines_of_file fn in
  let maybe_output_fn = CLI.get_string_opt ["-o"] args in
  let maybe_not_found_fn = CLI.get_string_opt ["-nf"] args in
  let not_found_names = ref StringSet.empty in
  let not_found =
    match maybe_not_found_fn with
    | None -> (fun name -> Log.warn "not found: %s" name)
    | Some _fn -> (fun name ->
        not_found_names := StringSet.add name !not_found_names
      ) in
  let force_db_creation = CLI.get_set_bool ["--force"] args in
  let no_index = CLI.get_set_bool ["--no-index"] args in
  let compress = CLI.get_set_bool ["-z"] args in
  let names_provider = match CLI.get_string_opt ["-names"] args with
    | Some names -> On_cli names
    | None ->
      let fn = CLI.get_string ["-f"] args in
      From_file fn in
  CLI.finalize ();
  let names = match names_provider with
    | On_cli names -> S.split_on_string names ~by:","
    | From_file fn -> LO.lines_of_file fn in
  let nb_names = L.length names in
  let out = match maybe_output_fn with
    | None -> stdout
    | Some output_fn ->
      if S.ends_with output_fn ".gz" then
        let cmd = sprintf "gzip -cq > %s" output_fn in
        Log.info "running: %s" cmd;
        Unix.open_process_out cmd
      else
        open_out_bin output_fn in
  (if no_index then
     begin
       let collected = match input_fns with
         | [input_fn] -> populate_ht names input_fn
         | [] -> failwith "Get_mol: --no-index requires -i"
         | _ :: _:: _ -> failwith "Get_mol: --no-index incompatible with -if" in
       let ht_len = Ht.length collected in
       if ht_len <> nb_names then
         Log.warn "found %d; expected %d" ht_len nb_names;
       L.iter (fun name ->
           try (* extract molecule *)
             let m = Ht.find collected name in
             fprintf out "%s" m
           with Not_found -> not_found name
         ) names
     end
   else
     begin
       let tmp_in_fn, tmp_out_fn, unlz4 =
         if compress then
           let tmp_in_fn = Fn.temp_file ~temp_dir:"/tmp" "get_mol_" ".lz4" in
           let tmp_out_fn = Fn.temp_file ~temp_dir:"/tmp" "get_mol_" "" in
           (tmp_in_fn, tmp_out_fn, lz4_string tmp_in_fn tmp_out_fn false)
         else
           ("", "", fun x -> x) in
       let dbs =
         L.map
           (db_open_or_create verbose force_db_creation compress)
           input_fns in
       L.iter (fun name ->
           try
             (* find containing db, if any *)
             let db = L.find (fun db -> DB.mem db name) dbs in
             (* extract molecule from it *)
             let m = unlz4 (DB.find db name) in
             fprintf out "%s" m
           with Not_found ->
             (* no db contains this molecule *)
             not_found name
         ) names;
       L.iter DB.close dbs;
       if tmp_in_fn <> "" then Unix.unlink tmp_in_fn;
       if tmp_out_fn <> "" then Unix.unlink tmp_out_fn
     end
  );
  (match maybe_output_fn with
   | Some _fn -> close_out out
   | None -> ()
  );
  (* dump not found names to file *)
  (match maybe_not_found_fn with
   | None -> ()
   | Some fn ->
     LO.with_out_file fn (fun out ->
         StringSet.iter (fun name ->
             fprintf out "%s\n" name
           ) !not_found_names
       )
  )

let () = main ()
