(* mini molecule module *)

open Printf

module A = BatArray
module IntSet = BatSet.Int
module L = MyList
module StringMap = BatMap.String

type t = { name: string;
           graph: Node.t array;
           matrix: int array array }

let get_name m = m.name

let get_graph m = m.graph

let nb_atoms m =
  A.length m.graph

let create name graph matrix =
  (* DEBUG *)
  (* let n = A.length matrix in
   * for i = 0 to n - 1 do
   *   for j = 0 to n - 1 do
   *     if j <> 0 then
   *       printf " %d" matrix.(i).(j)
   *     else
   *       printf "%d" matrix.(i).(j)
   *   done;
   *   printf "\n"
   * done; *)
  { name; graph; matrix }

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

let encode (max_height: int) (mol: t): (Atom_env.t * int) list =
  (* compute atom env. of given atom up to maximum height allowed *)
  let encode_atom (n_i: int): Atom_env.t =
    let rec loop height acc to_visit visited =
      if height > max_height || IntSet.is_empty to_visit then
        L.rev acc
      else
        let height' = height + 1 in
        let visited' = IntSet.union to_visit visited in
        let neighbors =
          IntSet.fold (fun x acc ->
              let neighbs = get_succs mol x in
              IntSet.union neighbs acc
            ) to_visit IntSet.empty in
        let to_visit' = IntSet.diff neighbors visited in
        let typs =
          IntSet.fold (fun x acc ->
              (get_typ mol x) :: acc
            ) to_visit' [] in
        (* canonicalize counted atom types *)
        let counted = Utls.list_uniq_count typs in
        let acc' = (height, counted) :: acc in
        loop height' acc' to_visit' visited'
    in
    let neighbors = get_succs mol n_i in
    let typ = get_typ mol n_i in
    (typ, loop 1 [] neighbors (IntSet.singleton n_i))
  in
  let nb_atoms = A.length mol.graph in
  let atom_indexes = L.range 0 `To (nb_atoms - 1) in
  (* canonicalize the encoding of the molecule by sorting its atom envs
     and counting duplicates *)
  let atom_envs = L.map encode_atom atom_indexes in
  Utls.list_uniq_count atom_envs
