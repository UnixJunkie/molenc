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
module Log = Dolog.Log
module Utls = Molenc.Utls

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_molecules a =
    create 1 Two_bands a
end

let verbose = ref false

(* FBR: several input files *)
(* FBR: input file so big that we cut it in chunks *)
(* FBR: integrate into FMGO *)
(* FBR: compress marshalled index? *)

(* FBR: in library module: index_many_from_files *)
(* FBR: in library module: nearest_in_many *)

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
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 50_000 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let i = ref 0 in
  let to_index, rest =
    let all_mols = FpMol.molecules_of_file input_fn in
    let l, r = L.takedrop csize all_mols in
    (ref l, ref r) in
  while !to_index <> [] do
    let chunk = A.of_list !to_index in
    let output_fn = sprintf "%s.%d.bst" input_fn !i in
    Log.info "creating %s" output_fn;
    let bst = Bstree.of_molecules chunk in
    Utls.save output_fn bst;
    let l, r = L.takedrop csize !rest in
    to_index := l;
    rest := r;
    incr i
  done

let () = main ()
