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

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let smiles_fn = Sys.argv.(1) in
  (* read each molecule *)
  LO.iter smiles_fn (fun line ->
      let smi, _name = S.split ~by:"\t" line in
      (* Log.info "smi: %s" smi; *)
      let mol_H =
        let mol = Rdkit.__init__ ~smi () in
        Rdkit.add_hydrogens mol () in
      let elements = Rdkit.get_elements mol_H () in
      (* get chemical formula *)
      let formula = formula_of_elements elements in
      Log.info "formula: %s" formula;
      (* encode to integer *)
      try
        let code = Formula.encode false formula in
        Log.info "code: %d" code;
        let formula' = Formula.decode code in
        Log.info "formula': %s" formula';
        assert(formula = formula')
      with Z.Overflow ->
        Log.error "code: overflow for %s" smi
    (* decode to chemical formula *)
    );
  ()

let () = main ()
