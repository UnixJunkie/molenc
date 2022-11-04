(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Work with DeepSMILES *)

open Printf

module A = Array
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module RNG = BatRandom.State
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

type mol = Rdkit.t

let mol_of_smiles (smi: string): mol =
  Rdkit.__init__ ~smi ()

let tok0 = "PAD"
let tok1 = "START"
let tok2 = "UNK"

let large_int = 1073741823 (* (2^30) - 1 *)

let lookup_in_dico dico word =
  try Ht.find dico word
  with Not_found ->
    let new_idx = Ht.length dico in
    Ht.add dico word new_idx;
    new_idx

let insert_reserved_tokens dico =
  let idx0 = lookup_in_dico dico tok0 in
  assert(idx0 = 0);
  let idx1 = lookup_in_dico dico tok1 in
  assert(idx1 = 1);
  let idx2 = lookup_in_dico dico tok2 in
  assert(idx2 = 3)

let main () =
  let start = Unix.gettimeofday () in  
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename>]: SMILES input file\n  \
              [-s <int>]: RNG seed\n  \
              [-n <int>]: DeepSMILES per input SMILES (default=1)\n  \
              -o <filename>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [--rand]: SMILES randomization (even if n=1)\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let maybe_seed = CLI.get_int_opt ["-s"] args in
  let n = CLI.get_int_def ["-n"] args 1 in
  let randomize = (CLI.get_set_bool ["--rand"] args) || n > 1 in
  let output_fn = CLI.get_string ["-o"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let rng = match maybe_seed with
    | Some s -> RNG.make [|s|] (* repeatable *)
    | None -> RNG.make_self_init () in
  let nb_molecules = ref 0 in
  let dummy = mol_of_smiles "C" in
  let dico = Ht.create 997 in (* enough to hold all DeepSMILES tokens *)
  insert_reserved_tokens dico;
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let line = input_line input in
          incr nb_molecules;
          let smi, _name = S.split line ~by:"\t" in
          let seed = RNG.int rng large_int in
          let dsmiles = Rdkit.get_deep_smiles dummy ~seed ~n ~randomize ~smi () in
          A.iter (Printf.fprintf output "%s\n") dsmiles
        done
      with End_of_file -> ()
    );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "encoding: %.1f molecule/s" ((float !nb_molecules) /. dt)

let () = main ()
