(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo, Japan.

   Counted atom pairs, with a new atom-typing scheme *)

open Printf

module A = Array
module Ht = BatHashtbl
module Rdkit = Molenc.Rdkit.Rdkit

let () = Py.initialize () (* because the Rdkit module uses Pyml *)

(* (atomic_num, nb_HA, nb_H, HA_used_val, formal_charge) cf. rdkit_wrapper.py *)
type atom = (int * int * int * int * int)
type pair = (atom * int * atom)

(* sorting the types for canonicalization *)
let create src dist dst =
  if compare src dst <= 0 then
    (src, dist, dst)
  else
    (dst, dist, src)

let string_of_atom (a,b,c,d,e) =
  sprintf "(%d,%d,%d,%d,%d)" a b c d e

let to_string (src, dist, dst) =
  Printf.sprintf "%s-%d-%s" (string_of_atom src) dist (string_of_atom dst)

let mol_of_smiles smi =
  Rdkit.__init__ ~smi ()

let type_atom mol i =
  let typ = Rdkit.type_atom mol ~i:i () in
  assert(A.length typ == 5);
  (typ.(0), typ.(1), typ.(2), typ.(3), typ.(4))

(* counted atom pairs *)
let encode_molecule mol =
  let n = Rdkit.get_num_atoms mol () in
  let ht = Ht.create (n * (n - 1) / 2) in
  for i = 0 to n - 1 do
    let src = type_atom mol i in
    for j = i to n - 1 do
      let dst = type_atom mol j in
      let dist = Rdkit.get_distance mol ~i ~j () in
      let feature = create src dist dst in
      let prev_count = Ht.find_default ht feature 0 in
      Ht.replace ht feature (prev_count + 1)
    done
  done;
  ht

let main () =
  let mol = mol_of_smiles "c1ccccc1" in
  let counted_pairs = encode_molecule mol in
  Ht.iter (fun k v -> printf "%s: %d\n" (to_string k) v) counted_pairs

let () = main ()
