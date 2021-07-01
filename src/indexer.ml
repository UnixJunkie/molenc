(* Copyright (C) 2021, Francois Berenger
   Tsuda laboratory, Tokyo university, Japan.

   Indexing of fingerprint encoded molecules into bisector trees
   stored on disk. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Utls = Molenc.Utls

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_molecules l =
    create 1 Two_bands (A.of_list l)
end

let verbose = ref false

(* FBR: several input files *)
(* FBR: integrate into FMGO *)

(* FBR: in library module: index_many_from_files *)
(* FBR: in library module: nearest_in_many *)

let read_one_chunk input_fn in_mol_count chunk_index csize input () =
  let res = ref [] in
  try
    for _i = 1 to csize do
      let line = input_line input in
      let mol = (!in_mol_count, line) in
      incr in_mol_count;
      res := mol :: !res
    done;
    let idx = !chunk_index in
    incr chunk_index;
    (idx, !res)
  with End_of_file ->
    if !res = [] then
      (Log.info "read %d from %s" !in_mol_count input_fn;
       raise Parany.End_of_input)
    else
      (* last chunk, maybe not full *)
      (!chunk_index, !res)

let index_one_chunk input_fn (i, chunk') =
  let chunk = L.rev_map (fun (i, line) -> FpMol.parse_one i line) chunk' in
  assert(i <= 9999);
  let output_fn = sprintf "%s.%04d.bst" input_fn i in
  Log.info "creating %s" output_fn;
  let bst = Bstree.of_molecules chunk in
  Utls.save output_fn bst;
  Utls.run_command (sprintf "gzip -f %s" output_fn)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: molecules input file\n  \
               -np <int>: nprocs (default=1)\n  \
               -c <int>: chunk size (molecules/bloc)\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 50_000 in
  if CLI.get_set_bool ["-v"] args then
    verbose := true;
  CLI.finalize (); (* ------------------------------------------------------ *)
  let chunk_count = ref 0 in
  let in_mol_count = ref 0 in
  Log.info "%d molecules in %s" (LO.count input_fn) input_fn;
  LO.with_in_file input_fn (fun input ->
      Parany.run nprocs
        ~demux:(read_one_chunk input_fn in_mol_count chunk_count csize input)
        ~work:(index_one_chunk input_fn)
        ~mux:(fun () -> ())
    )

let () = main ()
