(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* read files output by ./bin/type_atoms.py *)

module A = Array
module IntSet = BatSet.Int
module L = BatList

let read_one counter input =
  (* "#atoms:14 caffeine" *)
  let atoms_header = input_line input in
  let nb_atoms, mol_name =
    Scanf.sscanf atoms_header "#atoms:%d %s"
      (fun nb_atoms name -> (nb_atoms, name)) in
  (* read atoms *)
  let atom_lines = Utls.read_n_lines nb_atoms input in
  let atoms =
    L.map (fun l ->
        Scanf.sscanf l "%d %s"
          (fun _index typ -> PiEltHA.of_string typ)
      ) atom_lines in
  let atom_types = A.of_list atoms in
  (* read bonds header; like "#bonds:15" *)
  let bonds_header = input_line input in
  let nb_bonds = Scanf.sscanf bonds_header "#bonds:%d" (fun n -> n) in
  (* read bonds *)
  let bond_lines = Utls.read_n_lines nb_bonds input in
  let succs_table = A.make nb_atoms IntSet.empty in
  L.iter (fun bond_line ->
      Scanf.sscanf bond_line "%d %d" (fun start stop ->
          assert(start <> stop);
          (* we need to add the bond two times, because the molecular
             graph is undirected *)
          (* start -> stop *)
          succs_table.(start) <- IntSet.add stop succs_table.(start);
          (* stop -> start *)
          succs_table.(stop) <- IntSet.add start succs_table.(stop)
        )
    ) bond_lines;
  (* read distance matrix *)
  (* matrix header line *)
  let matrix_header = input_line input in
  let diameter = Scanf.sscanf matrix_header "#diameter:%d" (fun n -> n) in
  (* matrix' content *)
  let matrix_lines = Utls.read_n_lines nb_atoms input in
  let matrix = Array.make_matrix nb_atoms nb_atoms 0 in
  L.iteri (fun i line ->
      let dist_strings = BatString.split_on_string line ~by:" " in
      L.iteri (fun j str ->
          let d = int_of_string str in
          matrix.(i).(j) <- d
        ) dist_strings
    ) matrix_lines;
  let nodes =
    A.mapi (fun i typ ->
        Node.create typ succs_table.(i)
      ) atom_types in
  incr counter;
  Mini_mol.create mol_name nodes diameter matrix
