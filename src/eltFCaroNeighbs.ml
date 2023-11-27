(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Some new atom typing scheme for 2023 and later. *)

type t = { anum: int; (* atomic number *)
           fc: int; (* formal charge *)
           aro: bool; (* aromatic? *)
           neighbs: int array } (* count of bonded neighbors *)

(* FBR: to encode a molecule, rdkit must give us the topological distance matrix *)
(*      plus type each atom *)
