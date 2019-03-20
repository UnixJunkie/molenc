(* mini molecule module *)

module A = BatArray
module IntSet = BatSet.Int
module L = MyList
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

let to_string (m: t): string =
  let buff = Buffer.create 80 in
  Buffer.add_string buff m.name;
  Buffer.add_char buff '\n';
  A.iter (fun n ->
      Buffer.add_string buff (Node.to_string n);
      Buffer.add_char buff '\n';
    ) m.graph;
  Buffer.contents buff

let get_typ (m: t) (i: int) =
  Node.get_typ m.graph.(i)

let get_succs (m: t) (i: int) =
  Node.get_succs m.graph.(i)

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
