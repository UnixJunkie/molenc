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
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

type mol = Rdkit.t

let mol_of_smiles (smi: string): mol =
  Rdkit.__init__ ~smi ()

let read_one_line input () =
  try input_line input
  with End_of_file -> raise Parany.End_of_input

let process_one line =
  let smi, _name = S.split line ~by:"\t" in
  let mol = mol_of_smiles smi in
  failwith "not implemented yet"

let write_one total in_AD out = function
  | None -> incr total
  | Some line ->
    begin
      incr total;
      incr in_AD;
      fprintf out "%s\n" line
    end

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename>]: SMILES input file\n  \
              -o <filename>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let nb_molecules = ref 0 in
  let start = Unix.gettimeofday () in
  let atom_pairs = ref [] in
  LO.with_in_file input_fn (fun input ->
      failwith "not implemented yet"
      (* Parany.run nprocs ~demux:(read_one_line input) *)
      (*   ~work:process_one *)
      (*   ~mux:(record_one nb_molecules atom_pairs) *)
    );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "encoding: %.1f molecule/s" ((float !nb_molecules) /. dt)

let () = main ()
