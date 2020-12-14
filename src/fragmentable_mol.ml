
module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IS = BatSet.Int
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

open Printf

type atom = { pi_electrons: int;
              atomic_num: int;
              heavy_neighbors: int;
              formal_charge: int }

let dummy_atom = { pi_electrons = -1;
                   atomic_num = -1;
                   heavy_neighbors = -1;
                   formal_charge = -1 }

let atom_of_string s =
  Scanf.sscanf s "%d,%d,%d,%d" (fun i j k l ->
      { pi_electrons = i;
        atomic_num = j;
        heavy_neighbors = k;
        formal_charge = l }
    )

type bond_type = Single
               | Aromatic
               | Double
               | Triple

let bond_of_char = function
  | '-' -> Single
  | '~' -> Aromatic
  | '=' -> Double
  | '#' -> Triple
  | c -> failwith ("Fragmentable_mol.bond_of_char: unsupported: " ^
                   (String.make 1 c))

type bond = { start: int;
              btype: bond_type;
              stop: int }

let bond_of_string s =
  Scanf.sscanf s "%d %c %d" (fun start bchar stop ->
      { start;
        btype = bond_of_char bchar;
        stop }
    )

(* a molecule with fragmentation hints *)
type fragmentable =
  { atoms: atom array;
    bonds: bond array;
    cut_bonds: int array; (* indexes in the bonds array *)
    frag_hint: int } (* how many bonds are suggested to be broken *)

type attachment_point = { start: int; (* atom index *)
                          dest: atom } (* end-point allowed atom type *)

let dummy_attachment_point = { start = -1;
                               dest = dummy_atom }

type fragmented =
  { atoms: atom array;
    bonds: bond array;
    anchors: attachment_point array }

let parse_mol_header header =
  (* "^#atoms:19 NCGC00261763-01$" *)
  Scanf.sscanf header "#atoms:%d %s" (fun nb_atoms mol_name ->
      (nb_atoms, mol_name)
    )

let parse_bonds_header header =
  (* "^#bonds:33$" *)
  Scanf.sscanf header "#bonds:%d" (fun nb_bonds ->
      nb_bonds
    )

let parse_cut_bonds header =
  (* "^cut_bonds:%d:%d$" *)
  Scanf.sscanf header "#cut_bonds:%d:%d" (fun nb_cuttable cut_hint ->
      (nb_cuttable, cut_hint)
    )

let read_one_molecule (input: in_channel): fragmentable =
  let header = input_line input in
  let nb_atoms, mol_name = parse_mol_header header in
  Log.debug "%d %s" nb_atoms mol_name;
  (* read atoms *)
  let atoms =
    A.init nb_atoms (fun _i ->
        atom_of_string (input_line input)
      ) in
  (* read bonds *)
  let nb_bonds = parse_bonds_header (input_line input) in
  Log.debug "%d" nb_bonds;
  let bonds =
    A.init nb_bonds (fun _i ->
        bond_of_string (input_line input)
      ) in
  (* read cuttable_bonds *)
  let nb_cut_bonds, frag_hint = parse_cut_bonds (input_line input) in
  Log.debug "%d %d" nb_cut_bonds frag_hint;
  let cut_bonds =
    A.init nb_cut_bonds (fun _i ->
        int_of_string (input_line input)
      ) in
  (* return res *)
  { atoms;
    bonds;
    cut_bonds;
    frag_hint }

(* translate bonds to cut into attachment points *)
let attach_points (mol: fragmentable) (cut_bonds: int array)
  : attachment_point array =
  let n = A.length cut_bonds in
  let res = A.make (2*n) dummy_attachment_point in
  A.iteri (fun i to_cut ->
      let j = 2*i in
      let k = j + 1 in
      let bond = mol.bonds.(to_cut) in
      assert(bond.btype = Single); (* only those are "cuttable" *)
      res.(j) <- { start = bond.start;
                   dest = mol.atoms.(bond.stop) };
      res.(k) <- { start = bond.stop;
                   dest = mol.atoms.(bond.start) }
    ) cut_bonds;
  res

(* remove those bonds from the list of bonds;
   compute corresp. attachment points *)
let edit_bonds (mol: fragmentable) (to_cut: int array): fragmented =
  let anchors = attach_points mol to_cut in
  let prev_bonds = A.to_list mol.bonds in
  (* remove those bonds *)
  let curr_bonds =
    A.of_list (L.filteri (fun i _x -> not (A.mem i to_cut)) prev_bonds) in
  { atoms = mol.atoms;
    bonds = curr_bonds;
    anchors }

(* needed for computing connectivity later on *)
let directed_to_undirected_bonds (m: fragmented): bond array =
  let rev_bonds = A.map (fun orig ->
      { start = orig.stop; btype = orig.btype; stop = orig.start }
    ) m.bonds in
  A.append m.bonds rev_bonds

(* instant access to all successors of a node; query by its index *)
type succs_ht = (int, IS.t) Hashtbl.t

(* extract connectivity info *)
let compute_successors (m: fragmented): succs_ht =
  let res = Ht.create (A.length m.atoms) in
  let all_bonds = directed_to_undirected_bonds m in
  A.iter (fun (b: bond) ->
      let curr = b.start in
      let next = b.stop in
      try
        let prev_succs = Ht.find res curr in
        Ht.replace res curr (IS.add next prev_succs)
      with Not_found -> Ht.add res curr (IS.singleton next)
    ) all_bonds;
  res

(* extract whole connected fragment; starting from a single seed atom *)
let connected_component (seed: int) (succs: succs_ht): IS.t =
  let rec loop to_visit visited acc =
    if IS.is_empty to_visit then
      acc
    else
      let x, to_visit' = IS.pop to_visit in
      let visited' = IS.add x visited in
      let acc' = IS.add x acc in
      let successors = Ht.find succs x in
      let to_visit'' =
        IS.diff (IS.union to_visit' successors) visited' in
      let acc'' = IS.union acc' successors in
      loop to_visit'' visited' acc''
  in
  loop (IS.singleton seed) IS.empty IS.empty

(* cut a fragmented molecule into its disconnected fragments *)
let reconcile (_frag: fragmented): fragmented list =
  failwith "not implemented yet"

(* we fragment the molecule just once;
   larger molecules give rise to more fragments *)
let fragment_molecule rng m =
  let cuttable = A.copy m.cut_bonds in
  let nb_cuts = m.frag_hint in
  assert(A.length cuttable >= nb_cuts);
  A.shuffle ~state:rng cuttable;
  let to_cut = A.left cuttable nb_cuts in
  let _edited = edit_bonds m to_cut in
  (* this freshly cut molecule needs to be "reconciliated" *)
  (* compute connected components *)
  (* renumber nodes *)
  (* translate edges *)
  failwith "not implemented yet"

let main () =
  (* read in a file *)
  Log.(set_log_level DEBUG);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i molecules.txt\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  CLI.finalize();
  let _all_molecules =
    Utls.with_in_file input_fn (fun input ->
        let res, exn =
          L.unfold_exn (fun () -> read_one_molecule input) in
        (if exn <> End_of_file then raise exn);
        res
      ) in
  (* TODO: really fragment it then *)
  ()

let () = main ()
