(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * chemical formula integer encoder/decoder *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Formula = Molenc.Formula
module SMap = BatMap.String
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

let formula_of_elements (elts: string array): string =
  let elt2count =
    A.fold (fun acc elt ->
        let prev_count = SMap.find_default 0 elt acc in
        SMap.add elt (prev_count + 1) acc
      ) SMap.empty elts in
  let buff = Buffer.create 128 in
  SMap.iter (bprintf buff "%s%d") elt2count;
  Buffer.contents buff

(* FBR: run on whole ChEMBL to see how much whole molecules are encodable *)

type work_result = OK of string
                 | Overflow of string

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  (if argc = 1 then
     (eprintf "usage:\n  \
               %s -i in.smi\n  \
               -i <input.smi>: input molecules\n  \
               [-np <int>]: parallelize on NCORES (default=1)\n  \
               [-c <int>]: chunk size (default=50)\n"
        Sys.argv.(0);
      exit 1)
  );
  let smiles_fn = CLI.get_string ["-i"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 50 in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (* read each molecule *)
  LO.with_in_file smiles_fn (fun input ->
      Parany.run nprocs ~csize
        ~demux:(fun () ->
          try input_line input
          with End_of_file -> raise Parany.End_of_input
        )
        ~work:(fun line ->
          let smi, _name = S.split ~by:"\t" line in
          let mol_H =
            let mol = Rdkit.__init__ ~smi () in
            Rdkit.add_hydrogens mol () in
          let elements = Rdkit.get_elements mol_H () in
          (* get chemical formula *)
          let formula = formula_of_elements elements in
          (* Log.info "formula: %s" formula; *)
          (* encode to integer *)
          try
            let code = Formula.encode false formula in
            (* Log.info "code: %d" code; *)
            (* decode to chemical formula *)
            let formula' = Formula.decode code in
            (* Log.info "formula': %s" formula'; *)
            assert(formula = formula');
            OK line
          with Z.Overflow ->
            Overflow line
        )
        ~mux:(function
          | OK line -> printf "%s\n%!" line
          | Overflow line -> eprintf "%s\n%!" line
        )
    )

let () = main ()
