(* Copyright (C) 2021, Francois Berenger
   Tsuda laboratory, Tokyo university, Japan.

   Find name, Tanimoto-score and SMILES of nearest neighbor molecules. *)

open Printf

module A = BatArray
module Bstree = Molenc.Index.Bstree
module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Utls = Molenc.Utls

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  let default_block_size = ref 50_000 in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: molecules input file\n  \
               -ifs <filename>: file containing a list of files\n  \
               -np <int>: nprocs (default=1)\n  \
               -c <int>: chunk size (molecules/bloc; default=%d)\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0) !default_block_size;
      exit 1
    end;
  let input_fns =
    match (CLI.get_string_opt ["-i"] args,
           CLI.get_string_opt ["-ifs"] args) with
    | (None, None)
    | (Some _, Some _) -> failwith "provide either -i or -ifs"
    | (Some fn, None) -> [fn]
    | (None, Some fn) -> LO.lines_of_file fn in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args !default_block_size in
  if CLI.get_set_bool ["-v"] args then
    verbose := true;
  CLI.finalize (); (* ------------------------------------------------------ *)
  let chunk_count = ref 0 in
  let in_mol_count = ref 0 in
  L.iter (fun input_fn ->
      Log.info "%d molecules in %s" (LO.length input_fn) input_fn;
      LO.with_in_file input_fn (fun input ->
          Parany.run nprocs
            ~demux:(read_one_chunk input_fn in_mol_count chunk_count csize input)
            ~work:(index_one_chunk input_fn)
            ~mux:(fun () -> ())
        )
    ) input_fns

let () = main ()

(* MolIndex.nearest_neighbor_names_a ncores bst_fns to_annotate *)
  
