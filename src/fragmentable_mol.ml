
module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IS = BatSet.Int
module L = BatList
module Log = Dolog.Log
module RNG = BatRandom.State
module S = BatString
module Utls = Molenc.Utls

open Printf

type atom = { index: int;
              pi_electrons: int;
              atomic_num: int;
              heavy_neighbors: int;
              formal_charge: int }

type atom_type = int * int * int * int

type atom_type_pair = atom_type * atom_type

let dummy_atom_type = (-1, -1, -1, -1)

let get_atom_type at =
  (at.pi_electrons,
   at.atomic_num,
   at.heavy_neighbors,
   at.formal_charge)

let reindex_atom ht (at: atom): atom =
  { at with index = Ht.find ht at.index }

let compare_atom_indexes a1 a2 =
  compare a1.index a2.index

let dummy_atom = { index = -1;
                   pi_electrons = -1;
                   atomic_num = -1;
                   heavy_neighbors = -1;
                   formal_charge = -1 }

let atom_of_string s =
  Scanf.sscanf s "%d %d,%d,%d,%d"
    (fun a b c d e -> { index = a;
                        pi_electrons = b;
                        atomic_num = c;
                        heavy_neighbors = d;
                        formal_charge = e })

let string_of_atom a =
  sprintf "%d %d,%d,%d,%d"
    a.index
    a.pi_electrons
    a.atomic_num
    a.heavy_neighbors
    a.formal_charge

let string_of_atom_type (a, b, c, d) =
  sprintf "%d,%d,%d,%d" a b c d

type bond_type = Single
               | Aromatic
               | Double
               | Triple

let bond_of_char = function
  | '-' -> Single
  | ':' -> Aromatic
  | '=' -> Double
  | '#' -> Triple
  | c -> failwith ("Fragmentable_mol.bond_of_char: unsupported: " ^
                   (String.make 1 c))

let char_of_bond = function
  | Single -> '-'
  | Aromatic -> ':'
  | Double -> '='
  | Triple -> '#'

type bond = { start: int;
              btype: bond_type;
              stop: int }

let reindex_bond ht (b: bond): bond =
  { b with
    start = Ht.find ht b.start;
    stop = Ht.find ht b.stop }

let compare_bond_indexes b1 b2 =
  compare (b1.start, b1.stop) (b2.start, b2.stop)

let bond_of_string s =
  Scanf.sscanf s "%d %c %d"
    (fun start bchar stop -> { start;
                               btype = bond_of_char bchar;
                               stop })

let string_of_bond b =
  sprintf "%d %c %d" b.start (char_of_bond b.btype) b.stop

(* a molecule with fragmentation hints *)
type fragmentable =
  { name: string;
    atoms: atom array;
    bonds: bond array;
    cut_bonds: int array; (* indexes in the bonds array *)
    frag_hint: int } (* how many bonds are suggested to be broken *)

type molecule =
  { name: string;
    atoms: atom array;
    bonds: bond array }

let write_one_molecule out mol =
  fprintf out "#atoms:%d %s\n" (A.length mol.atoms) mol.name;
  A.iter (fun a -> fprintf out "%s\n" (string_of_atom a)) mol.atoms;
  fprintf out "#bonds:%d\n" (A.length mol.bonds);
  A.iter (fun b -> fprintf out "%s\n" (string_of_bond b)) mol.bonds

type anchor = { src_typ: atom_type; (* src atom type *)
                start: int; (* src atom index *)
                dst_typ: atom_type } (* end-point allowed atom type *)

let get_anchor_type a =
  (a.src_typ, a.dst_typ)

let get_flipped_anchor_type a =
  (a.dst_typ, a.src_typ)

let flip_pair (a, b) =
  (b, a)

let compatible_anchors a1 a2 =
  (a1.src_typ = a2.dst_typ) && (a1.dst_typ = a2.src_typ)

let reindex_anchor ht (a: anchor): anchor =
  { a with start = Ht.find ht a.start }

let get_anchor_start (a: anchor): int =
  a.start

let dummy_anchor = { src_typ = dummy_atom_type;
                     start = -1;
                     dst_typ = dummy_atom_type }

let string_of_anchor a =
  sprintf "%s %d %s"
    (string_of_atom_type a.src_typ)
    a.start
    (string_of_atom_type a.dst_typ)

let anchor_of_string s =
  Scanf.sscanf s "%d,%d,%d,%d %d %d,%d,%d,%d"
    (fun a b c d start i j k l ->
       { src_typ = (a, b, c, d);
         start;
         dst_typ = (i, j, k, l) })

type fragment =
  { name: string;
    atoms: atom array;
    bonds: bond array;
    anchors: anchor array }

(* remove one compatible (randomly selected) anchor from
   the fragments' anchors *)
let rand_remove_anchor rng (src_typ, dst_typ) frag =
  let fake_anchor = { src_typ; start = -1; dst_typ } in
  let indexes =
    A.fold_righti (fun i x acc ->
        if compatible_anchors fake_anchor x then
          i :: acc
        else
          acc
      ) frag.anchors [] in
  match indexes with
  | [] -> failwith "Fragmentable_mol.rand_remove_anchor: should never happen"
  | _ ->
    (* draw one *)
    let i = Utls.array_random_elt rng (A.of_list indexes) in
    let anchors' = Utls.array_without_elt_at i frag.anchors in
    { frag with anchors = anchors' }

let rand_select_compat_anchor rng (src_typ, dst_typ) frag =
  let fake_anchor = { src_typ; start = -1; dst_typ } in
  let indexes =
    A.fold_righti (fun i x acc ->
        if compatible_anchors fake_anchor x then
          i :: acc
        else
          acc
      ) frag.anchors [] in
  match indexes with
  | [] -> assert(false)
  | _ -> Utls.array_random_elt rng (A.of_list indexes)

(* tree data structure to construct a molecule *)
type mol_tree = Branch of fragment * mol_tree list
              | Leaf of fragment

let get_name = function
  | Branch (b, _) -> b.name
  | Leaf l -> l.name

exception No_such_anchor of string

(* if there are several adequate anchors, we need to choose
 * one randomly (and allow for repeatibility) so that the fragment
 * linking process is not biased *)
let rand_select_anchor
    rng (a: anchor) (tree: mol_tree): int * mol_tree =
  let frag = match tree with
    | Leaf l -> l
    | Branch (core, _branches) -> core in
  (* find anchor indexes with the right type *)
  let indexes =
    A.fold_righti (fun i x acc ->
        if compatible_anchors a x then
          i :: acc
        else
          acc
      ) frag.anchors [] in
  match indexes with
  | [] -> raise (No_such_anchor
                   (sprintf "no compatible anchor for (%s) in %s"
                      (string_of_anchor a) frag.name))
  | _ ->
    (* draw one *)
    let i = Utls.array_random_elt rng (A.of_list indexes) in
    let target = frag.anchors.(i) in
    let anchors' = Utls.array_without_elt_at i frag.anchors in
    (* return target atom index and update anchors list *)
    let frag' = { frag with anchors = anchors' } in
    let tree' = match tree with
      | Leaf _ -> Leaf frag'
      | Branch (_core, branches) -> Branch (frag', branches) in
    (target.start, tree')

let append_index_to_frag_name (index: int) (frag: fragment): fragment =
  { frag with name = sprintf "%s_%03d" frag.name index }

(* renumber all atoms, bonds and anchors *)
let reindex offset frag =
  let n = A.length frag.atoms in
  let ht = Ht.create n in
  A.iteri (fun i at ->
      let old_index = at.index in
      let new_index = i + !offset in
      Ht.add ht old_index new_index
    ) frag.atoms;
  let atoms = A.map (reindex_atom ht) frag.atoms in
  let bonds = A.map (reindex_bond ht) frag.bonds in
  let anchors = A.map (reindex_anchor ht) frag.anchors in
  offset := !offset + n; (* update offset *)
  { name = frag.name; atoms; bonds; anchors }

let get_anchor_types (f: fragment): atom_type_pair list =
  (* in the right order *)
  A.fold_right (fun x acc ->
      (get_anchor_type x) :: acc
    ) f.anchors []

let get_anchor_starts (f: fragment): int list =
  (* in the right order *)
  A.fold_right (fun x acc ->
      (get_anchor_start x) :: acc
    ) f.anchors []

let write_one_fragment out name index frag =
  fprintf out "#atoms:%d %s_f%02d\n" (A.length frag.atoms) name !index;
  incr index;
  A.iter (fun a -> fprintf out "%s\n" (string_of_atom a)) frag.atoms;
  fprintf out "#bonds:%d\n" (A.length frag.bonds);
  A.iter (fun b -> fprintf out "%s\n" (string_of_bond b)) frag.bonds;
  fprintf out "#anchors:%d\n" (A.length frag.anchors);
  A.iter (fun a -> fprintf out "%s\n" (string_of_anchor a)) frag.anchors

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
  (* "^#cut_bonds:%d:%d$" *)
  Scanf.sscanf header "#cut_bonds:%d:%d" (fun nb_cuttable cut_hint ->
      (nb_cuttable, cut_hint)
    )

let parse_anchors header =
  (* "^#anchors:2$" *)
  Scanf.sscanf header "#anchors:%d" (fun x -> x)

let read_one_molecule (input: in_channel): fragmentable =
  let header = input_line input in
  let nb_atoms, name = parse_mol_header header in
  Log.debug "%d %s" nb_atoms name;
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
  { name;
    atoms;
    bonds;
    cut_bonds;
    frag_hint }

(* ---
#atoms:6 NCGC00261552-01_f02
10 1,6,3,0
11 1,6,2,0
12 1,6,2,0
13 1,6,2,0
14 1,6,2,0
15 1,6,2,0
#bonds:6
10 ~ 11
11 ~ 12
12 ~ 13
13 ~ 14
14 ~ 15
15 ~ 10
#anchors:1
10 3 0,6,3,0
--- *)
let read_one_fragment (input: in_channel): fragment =
  (* read atoms *)
  let nb_atoms, name = parse_mol_header (input_line input) in
  let atoms =
    A.init nb_atoms (fun _i ->
        atom_of_string (input_line input)
      ) in
  (* read bonds *)
  let bonds =
    let nb_bonds = parse_bonds_header (input_line input) in
    A.init nb_bonds (fun _i ->
        bond_of_string (input_line input)
      ) in
  (* read attachment points *)
  let anchors =
    let nb_anchors = parse_anchors (input_line input) in
    A.init nb_anchors (fun _i ->
        anchor_of_string (input_line input)
      ) in
  (* return res *)
  Log.debug "%d %d %d" (A.length atoms) (A.length bonds) (A.length anchors);
  { name; atoms; bonds; anchors }

(* translate bonds to cut into attachment points *)
let anchor_points (mol: fragmentable) (cut_bonds: int array): anchor array =
  let n = A.length cut_bonds in
  let res = A.make (2*n) dummy_anchor in
  A.iteri (fun i to_cut ->
      let j = 2*i in
      let k = j + 1 in
      let bond = mol.bonds.(to_cut) in
      assert(bond.btype = Single); (* only those are "cuttable" *)
      res.(j) <- { src_typ = get_atom_type mol.atoms.(bond.start);
                   start = bond.start;
                   dst_typ = get_atom_type mol.atoms.(bond.stop) };
      res.(k) <- { src_typ = get_atom_type mol.atoms.(bond.stop);
                   start = bond.stop;
                   dst_typ = get_atom_type mol.atoms.(bond.start) }
    ) cut_bonds;
  res

(* remove those bonds from the list of bonds;
   compute corresp. attachment points *)
let edit_bonds (mol: fragmentable) (to_cut: int array): fragment =
  let anchors = anchor_points mol to_cut in
  let prev_bonds = A.to_list mol.bonds in
  (* remove those bonds *)
  let curr_bonds =
    A.of_list (L.filteri (fun i _x -> not (A.mem i to_cut)) prev_bonds) in
  (* inherit parent molecule name; an index will be appended later *)
  { name = mol.name;
    atoms = mol.atoms;
    bonds = curr_bonds;
    anchors }

(* needed for computing connectivity later on *)
let directed_to_undirected_bonds (m: fragment): bond array =
  let rev_bonds = A.map (fun orig ->
      { start = orig.stop; btype = orig.btype; stop = orig.start }
    ) m.bonds in
  A.append m.bonds rev_bonds

(* instant access to all successors of a node; query by its index *)
type succs_ht = (int, IS.t) Hashtbl.t

(* extract connectivity info *)
let compute_successors (m: fragment): succs_ht =
  let res = Ht.create (A.length m.atoms) in
  let all_bonds = directed_to_undirected_bonds m in
  A.iter (fun (b: bond) ->
      let curr = b.start in
      let next = b.stop in
      let prev_succs = Ht.find_default res curr (IS.singleton next) in
      Ht.replace res curr (IS.add next prev_succs)
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
      (* isolated heavy atom (smallest possible fragment) --> IS.empty *)
      let successors = Ht.find_default succs x IS.empty in
      let to_visit'' =
        IS.diff (IS.union to_visit' successors) visited' in
      let acc'' = IS.union acc' successors in
      loop to_visit'' visited' acc''
  in
  loop (IS.singleton seed) IS.empty IS.empty

let enforce_conn_comp (comp: IS.t) (frag: fragment): fragment =
  { name = frag.name;
    atoms = A.filteri (fun i _x -> IS.mem i comp) frag.atoms;
    bonds =
      A.filter (fun (b: bond) ->
          (IS.mem b.start comp) && (IS.mem b.stop comp)
        ) frag.bonds;
    anchors = A.filter (fun a -> IS.mem a.start comp) frag.anchors }

(* cut fragmented molecule into its fragments *)
let reconcile (frag: fragment): fragment list =
  let succs = compute_successors frag in
  let processed = ref IS.empty in
  let count = ref 0 in
  A.fold_left (fun acc (a: anchor) ->
      if IS.mem a.start !processed then
        acc
      else
        let conn = connected_component a.start succs in
        let res = enforce_conn_comp conn frag in
        let res' = append_index_to_frag_name !count res in
        incr count;
        processed := IS.union conn !processed;
        res' :: acc
    ) [] frag.anchors

(* we fragment the molecule just once;
   larger molecules give rise to more fragments *)
let fragment_molecule rng m =
  let cuttable = A.copy m.cut_bonds in
  let nb_cuts = m.frag_hint in
  assert(A.length cuttable >= nb_cuts);
  A.shuffle ~state:rng cuttable;
  let to_cut = A.left cuttable nb_cuts in
  let edited = edit_bonds m to_cut in
  (* this freshly cut molecule needs to be "reconciliated" *)
  reconcile edited

type mode = Fragment of string * string (* (in_mols_fn, out_frags_fn) *)
          | Assemble of string * string (* (in_frags_fn, out_mols_fn) *)

(* organize fragments in such a way that drawing some is efficient *)
let organize_fragments frags_a =
  let ht = Ht.create (A.length frags_a) in
  (* list all atom types of attachment points *)
  (* for each atom type pair, count the number of fragments *)
  A.iter (fun frag ->
      A.iter (fun a ->
          let k = get_flipped_anchor_type a in
          let prev_frags = Ht.find_default ht k [] in
          Ht.replace ht k (frag :: prev_frags)
        ) frag.anchors
    ) frags_a;
  (* create Ht: atom_type pair -> fragments array *)
  Ht.map (fun _typs frags_lst ->
      A.of_list frags_lst
    ) ht

(* prepare for final molecule *)
let reindex_mol_tree t =
  let count = ref 0 in
  let rec loop = function
    | Leaf leaf -> Leaf (reindex count leaf)
    | Branch (core, branches) -> Branch (reindex count core,
                                         L.map loop branches) in
  loop t

(* extract all fragment names from the tree *)
let string_of_mol_tree t =
  let buff = Buffer.create 1024 in
  let rec loop = function
    | Leaf l -> bprintf buff "%s\n" l.name
    | Branch (core, branches) ->
      let branch_names = L.map get_name branches in
      L.iter (fun name ->
          bprintf buff "%s -> %s\n" core.name name
        ) branch_names;
      L.iter loop branches in
  loop t;
  Buffer.contents buff

let molecule_from_tree rng (name_prfx: string) (t: mol_tree): molecule =
  let rec loop (atoms, bonds, names) = function
    | Leaf leaf ->
      (* alredy linked to its parent fragment *)
      (Utls.prepend_list_with_array leaf.atoms atoms,
       Utls.prepend_list_with_array leaf.bonds bonds,
       leaf.name :: names)
    | Branch (core, branches) ->
      let atoms' = Utls.prepend_list_with_array core.atoms atoms in
      let bonds' = Utls.prepend_list_with_array core.bonds bonds in
      let names' = core.name :: names in
      let starts = get_anchor_starts core in
      let anchors = A.to_list core.anchors in
      let target_indexes_trees =
        L.map2 (rand_select_anchor rng) anchors branches in
      let stops, branches' = L.split target_indexes_trees in
      (* add new bonds *)
      let new_bonds =
        L.rev_map2 (fun start stop ->
            { start; btype = Single; stop }
          ) starts stops in
      let bonds'' = L.rev_append new_bonds bonds' in
      let acc = (atoms', bonds'', names') in
      L.fold_left loop acc branches' in
  let atoms_acc, bonds_acc, names_acc =
    try loop ([], [], []) t with
    | No_such_anchor s ->
      begin
        Log.fatal "%s" s;
        let tree_str = string_of_mol_tree t in
        Log.fatal "\n%s" tree_str;
        exit 1
      end in
  (* sort atoms, bonds and names to ease verification *)
  let atoms = A.of_list atoms_acc in
  let bonds = A.of_list bonds_acc in
  let names = L.sort BatString.compare names_acc in
  let frag_names = S.concat "," names in
  let name = sprintf "%s_%s" name_prfx frag_names in
  A.sort compare_atom_indexes atoms;
  A.sort compare_bond_indexes bonds;
  { name; atoms; bonds }

(* (* draw all necessary fragments *)
 * (* WARNING: not tail rec *)
 * let rec get_frag_for_anchor_type rng state typ frags_ht =
 *   let candidate = Utls.array_random_elt rng (Ht.find frags_ht typ) in
 *   if state = Grow && A.length candidate.anchors = 1 then
 *     Leaf candidate (* terminal fragment *)
 *   else
 *     let compat_typ = flip_pair typ in
 *     (* extract types of anchor points *)
 *     let all_types = get_anchor_types candidate in
 *     (* remove the one we just connected, if applicable *)
 *     let candidate', rem_types = match state with
 *       | Seed -> (candidate, all_types)
 *       | Grow -> (rand_remove_anchor rng typ candidate,
 *                  Utls.list_really_remove_one all_types compat_typ) in
 *     (* rec call *)
 *     Branch (candidate',
 *             L.map (fun typ' ->
 *                 get_frag_for_anchor_type rng Grow typ' frags_ht
 *               ) rem_types) *)

(* let connect_fragments rng name_prfx frags_a frags_ht =
 *   (* draw uniformly seed fragment *)
 *   let seed_frag = Utls.array_random_elt rng frags_a in
 *   let anchor = Utls.array_random_elt rng seed_frag.anchors in
 *   let typ = get_anchor_type anchor in
 *   let fragments = get_frag_for_anchor_type rng Seed typ frags_ht in
 *   let reindexed = reindex_mol_tree fragments in
 *   molecule_from_tree rng name_prfx reindexed *)

type state = Seed of fragment (* fragment with no ancestor in the tree *)
           | Grow of int * fragment (* fragment with already one anchor used *)

(* draw fragments to build molecule starting from seed fragment *)
let grow rng seed_frag frags_ht =
  let rec loop = function
    | Seed frag ->
      let anchors = A.to_list frag.anchors in
      let children_frags =
        L.map (fun anchor ->
            let typ = get_anchor_type anchor in
            (* draw one compatible fragment *)
            let cand_frag = Utls.array_random_elt rng (Ht.find frags_ht typ) in
            (* random choose index of compatible anchor in cand_frag *)
            let i = rand_select_compat_anchor rng typ cand_frag in
            Grow (i, cand_frag)
          ) anchors in
      Branch (frag, L.map loop children_frags)
    | Grow (i, frag) ->
      let anchors = frag.anchors in
      let n = A.length anchors in
      assert(n > 0);
      if n = 1 then
        Leaf frag
      else
        let anchors' = A.to_list (Utls.array_without_elt_at i anchors) in
        let children_frags =
          L.map (fun anchor ->
              let typ = get_anchor_type anchor in
              (* draw one compatible fragment *)
              let cand_frag = Utls.array_random_elt rng (Ht.find frags_ht typ) in
              (* random choose index of compatible anchor in cand_frag *)
              let j = rand_select_compat_anchor rng typ cand_frag in
              Grow (j, cand_frag)
            ) anchors' in
        Branch (frag, L.map loop children_frags) in
  loop (Seed seed_frag)
    
let connect_fragments rng name_prfx frags_a frags_ht =
  (* draw uniformly seed fragment *)
  let seed_frag = Utls.array_random_elt rng frags_a in
  let fragments = grow rng seed_frag frags_ht in
  let reindexed = reindex_mol_tree fragments in
  molecule_from_tree rng name_prfx reindexed

(* FBR: check some generated mols; create regression test suite *)

let main () =
  Log.(set_log_level INFO);
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-im <molecules.to_frag>]: input file; molecules w/ typed atoms\n  \
              and fragmenting hints (generated by molenc_frag.py)\n  \
              [-of <molecules.frags>]: computed fragments output file\n \
              [-if <molecules.frags>]: input fragments file\n  \
              [-om <molecules.txt>]: generated molecules output file\n  \
              [-s <rng_seed:int>]: for reproducibility\n  \
              [-p <int>]: number of fragmentation passes (default=1)\n  \
              [-n <int>]: nb. molecules to generate (default=1)\n  \
              [-v]: verbose log\n"
       Sys.argv.(0);
     exit 1);
  let maybe_in_mols_fn = CLI.get_string_opt ["-im"] args in
  let maybe_out_frags_fn = CLI.get_string_opt ["-of"] args in
  let maybe_in_frags_fn = CLI.get_string_opt ["-if"] args in
  let maybe_out_mols_fn = CLI.get_string_opt ["-om"] args in
  let nb_frag_passes = CLI.get_int_def ["-p"] args 1 in
  let nb_mols = CLI.get_int_def ["-n"] args 1 in
  let rng = match CLI.get_int_opt ["-s"] args with
    | Some s -> RNG.make [|s|] (* repeatable *)
    | None -> RNG.make_self_init () in
  if CLI.get_set_bool ["-v"] args then
    Log.(set_log_level DEBUG)
  ;
  CLI.finalize();
  let mode = match maybe_in_mols_fn, maybe_out_frags_fn with
    | (Some ifn, Some ofn) -> Fragment (ifn, ofn)
    | (None, None) ->
      begin match maybe_in_frags_fn, maybe_out_mols_fn with
        | (Some ifn, Some ofn) -> Assemble (ifn, ofn)
        | _ -> failwith "provide either -im and -of XOR -if and -om"
      end
    | _ -> failwith "provide either -im and -of XOR -if and -om" in
  match mode with
  | Fragment (input_fn, output_fn) ->
    (* read in molecules with fragment hints *)
    let all_molecules =
      Utls.with_in_file input_fn (fun input ->
          let res, exn = L.unfold_exn (fun () -> read_one_molecule input) in
          (if exn <> End_of_file then raise exn);
          res
        ) in
    (* fragment them *)
    (* parallelizable if needed; but already super fast; 6567 molecule/s *)
    Utls.with_out_file output_fn (fun out ->
        L.iter (fun mol ->
            let index = ref 0 in
            for _i = 1 to nb_frag_passes do
              let frags = fragment_molecule rng mol in
              L.iter (write_one_fragment out mol.name index) frags
            done
          ) all_molecules
      )
  | Assemble (input_fn, output_fn) ->
    (* read in fragments *)
    let dt0, all_fragments =
      Utls.time_it (fun () ->
          Utls.with_in_file input_fn (fun input ->
              let res, exn = L.unfold_exn (fun () -> read_one_fragment input) in
              (if exn <> End_of_file then raise exn);
              A.of_list res
            )
        ) in
    let nb_fragments = A.length all_fragments in
    Log.info "read %d fragments in %1.2fs" nb_fragments dt0;
    let dt1, frags_ht =
      Utls.time_it (fun () -> organize_fragments all_fragments) in
    Log.info "compiled fragments in %1.2fs" dt1;
    Utls.with_out_file output_fn (fun out ->
        let dt2, () =
          Utls.time_it (fun () ->
              for i = 1 to nb_mols do
                let name_prfx = sprintf "genmol_%06d" i in
                let mol = connect_fragments rng name_prfx all_fragments frags_ht in
                write_one_molecule out mol
              done
            ) in
        Log.info "gen %d mols in %1.2fs" nb_mols dt2
      )

let () = main ()
