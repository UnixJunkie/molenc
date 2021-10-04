(* Copyright (C) 2021, Francois Berenger
   Tsuda laboratory,
   Graduate School of Frontier Sciences,
   The University of Tokyo,
   5-1-5 Kashiwa-no-ha,
   Kashiwa, Chiba 277-8561, Japan.

   Find name, Tanimoto-score and SMILES of nearest neighbor molecules.
   Output a valid SMILES file. *)

open Printf

module A = BatArray
module Bstree = Molenc.Index.Bstree
module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module MolIndex = Molenc.Index
module S = BatString
module Utls = Molenc.Utls

let split_smiles_line l =
  (* Expect '\t' separated SMILES *)
  S.split l ~by:"\t"

let process_smiles_file name2smi fn =
  LO.iteri fn (fun i line ->
      let smi, name = split_smiles_line line in
      (if Ht.mem name2smi name then
         let () = Log.fatal "Finder.process_smiles_file: \
                             already seen molecule name: %s" name in
         exit 1
       else
         Ht.add name2smi name smi
      );
      if (i mod 1000) = 0 then
        printf "Loaded molecules: %d\r%!" (Ht.length name2smi)
    )

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename.AP>: encoded molecules input file\n  \
               --bst-fns <fn1[,fn2[,fn3...]]>: list of BST index files\n  \
               --smi-fns <fn1[,fn2[,fn3...]]>: list of SMILES files\n  \
               -np <int>: nprocs (default=1)\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let bst_fns = S.split_on_char ',' (CLI.get_string ["--bst-fns"] args) in
  let smi_fns = S.split_on_char ',' (CLI.get_string ["--smi-fns"] args) in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (* populate the name to SMILES LUT *)
  let name2smi = Ht.create 1_000_000 in
  L.iter (process_smiles_file name2smi) smi_fns;
  let encoded_molecules_in =
    A.of_list (Molenc.FpMol.molecules_of_file input_fn) in
  LO.with_out_file output_fn (fun out ->
      let fp_name_dists =
        MolIndex.nearest_neighbor_names_a
          nprocs bst_fns encoded_molecules_in in
      A.iter (fun (_fp, name, dist) ->
          let smi = Ht.find name2smi name in
          let tani = 1.0 -. dist in
          fprintf out "%s\t%s_%.2f\n" smi name tani
        ) fp_name_dists
    )
  
let () = main ()
