(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* mini molecule module *)

module A = BatArray
module IntSet = BatSet.Int
module Ht = BatHashtbl
module L = BatList
module StringMap = BatMap.String

type t = { name: string;
           graph: Node.t array;
           diameter: int;
           matrix: int array array }

let get_name m = m.name

let get_graph m = m.graph

let nb_atoms m =
  A.length m.graph

let create name graph diameter matrix =
  { name; graph; diameter; matrix }

let get_typ (m: t) (i: int) =
  Node.get_typ (A.unsafe_get m.graph i)

let get_succs (m: t) (i: int) =
  Node.get_succs (A.unsafe_get m.graph i)

(* list (sorted-uniq-counted) atom types of all atoms
   at given distance from center atom *)
let types_at_distance (center: int) (curr_height: int) (mol: t) =
  let matrix_line = mol.matrix.(center) in
  let unsorted =
    A.fold_lefti (fun acc i x ->
        if x = curr_height then
          (get_typ mol i) :: acc
        else
          acc
      ) [] matrix_line in
  (* layer at 'curr_height' *)
  (curr_height, Utls.list_uniq_count unsorted)

let encode (max_height: int) (mol: t): (Atom_env.t * int) list =
  (* compute atom envs. of given atom up to maximum height allowed *)
  (* we cannot go deeper than 'maxi' on this molecule *)
  let maxi = min max_height mol.diameter in
  let encode_atom (n_i: int): Atom_env.t =
    let depths = L.range 0 `To maxi in
    let layers =
      L.map (fun height ->
          types_at_distance n_i height mol
        ) depths in
    (* non empty layers *)
    L.filter (fun (_h, typs) -> typs <> []) layers
  in
  let nb_atoms = A.length mol.graph in
  let atom_indexes = L.range 0 `To (nb_atoms - 1) in
  (* canonicalize the encoding of the molecule by sorting its atom envs
     and counting duplicates *)
  let atom_envs = L.map encode_atom atom_indexes in
  Utls.list_uniq_count atom_envs

(* encode the molecule to counted atom pairs *)
let atom_pairs (mol: t): (Atom_pair.t * int) list =
  let n = nb_atoms mol in
  assert(n >= 1); (* at least one heavy atom *)
  let max_nb_pairs = max 1 (n * (n - 1) / 2) in
  let pair2count = Ht.create max_nb_pairs in
  for i = 0 to n - 1 do
    let type_i = get_typ mol i in
    for j = i to n - 1 do
      let type_j = get_typ mol j in
      let dist = A.unsafe_get (A.unsafe_get mol.matrix i) j in
      let pair = Atom_pair.create type_i type_j dist in
      let prev_count = Ht.find_default pair2count pair 0 in
      Ht.replace pair2count pair (prev_count + 1)
    done;
  done;
  let pair_counts = Ht.bindings pair2count in
  L.sort compare pair_counts (* canonicalization by sorting *)
