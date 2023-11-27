(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * Some new atom typing scheme for 2023 and later. *)

type t = { anum: int; (* atomic number *)
           fc: int; (* formal charge *)
           aro: bool; (* aromatic? *)
           neighbs: char array } (* count bonded neighbors; maximum count is 6 (max. valence of sulfur) *)

(* FBR: to encode a molecule, rdkit must give us the topological distance matrix *)
(*      plus type each atom *)

let index_of_anum = function
  | 6  (* C  *) -> 1
  | 1  (* H  *) -> 2
  | 7  (* N  *) -> 3
  | 8  (* O  *) -> 4
  | 15 (* P  *) -> 5
  | 16 (* S  *) -> 6
  | 9  (* F  *) -> 7
  | 17 (* Cl *) -> 8
  | 35 (* Br *) -> 9
  | 53 (* I  *) -> 10
  | x ->
    let () = Log.warn "index_of_anum: unsupported anum: %d" x in
    0 (* future proof; don't crash at least *)

let char_of_count = function
  | 0 -> 'a'
  | 1 -> 'b'
  | 2 -> 'c'
  | 3 -> 'd'
  | 4 -> 'e'
  | 5 -> 'f'
  | x ->
    let () = Log.fatal "char_of_count: unsupported count: %d" x in
    assert(false)
