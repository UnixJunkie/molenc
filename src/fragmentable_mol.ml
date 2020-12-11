
module A = BatArray
module CLI = Minicli.CLI
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

open Printf
               
type atom = { pi_electrons: int;
              atomic_num: int;
              heavy_neighbors: int;
              formal_charge: int }

(* unused *)
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

let read_one_molecule input =
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
