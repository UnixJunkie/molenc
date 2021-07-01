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
(* FBR: input file so big that we cut it in chunks *)
(* FBR: integrate into FMGO *)
(* FBR: compress marshalled index? *)

(* FBR: in library module: index_many_from_files *)
(* FBR: in library module: nearest_in_many *)

let read_one_chunk in_count csize input =
  let res = ref [] in
  for _i = 1 to csize do
    let line = input_line input in
    let mol = FpMol.parse_one !in_count line in
    incr in_count;
    res := mol :: !res
  done;
  !res

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
  let in_count = ref 0 in
  LO.with_in_file input_fn (fun input ->
      try
        while true do
          let output_fn = sprintf "%s.%d.bst" input_fn !i in
          incr i;
          Log.info "creating %s" output_fn;
          let chunk = read_one_chunk in_count csize input in
          let bst = Bstree.of_molecules chunk in
          Utls.save output_fn bst
        done
      with End_of_file ->
        Log.info "From %s, read %d" input_fn !in_count
    )

let () = main ()
