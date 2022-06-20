(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Read molecule name, element symbols and their 3D coordinates from a .sdf file
   holding 3D conformers. *)

module A = BatArray
module L = BatList
module Log = Dolog.Log
module S = BatString
module V3 = Vector3

type atoms_3D = { name: string;
                  elements: int array;
                  coords: Vector3.t array;
                  (* just which atom is connected to which other;
                     no bond order info *)
                  bonds: int list array }

type encoded_atom = float array array (* shape: (nb_dx, nb_channels) *)

let anum_of_symbol = function
  | "C" -> 6
  | "H" -> 1
  | "N" -> 7
  | "O" -> 8
  | "P" -> 15
  | "S" -> 16
  | "F" -> 9
  | "Cl" -> 17
  | "Br" -> 35
  | "I" -> 53
  | unk ->
    let () = Log.warn "unsuported elt: %s" unk in
    -1

let symbol_of_anum = function
 |  6  -> "C"
 |  1  -> "H"
 |  7  -> "N"
 |  8  -> "O"
 | 15  -> "P"
 | 16  -> "S"
 |  9  -> "F"
 | 17  -> "Cl"
 | 35  -> "Br"
 | 53  -> "I"
 | -1  -> "_" (* unsupported elt. *)
 | _ -> assert(false)

let nb_channels = 10

let channel_of_anum = function
 |  6  -> 0 (* "C" *)
 |  1  -> 1 (* "H" *)
 |  7  -> 2 (* "N" *)
 |  8  -> 3 (* "O" *)
 | 15  -> 4 (* "P" *)
 | 16  -> 5 (* "S" *)
 |  9  -> 6 (* "F" *)
 | 17  -> 7 (* "Cl" *)
 | 35  -> 8 (* "Br" *)
 | 53  -> 9 (* "I" *)
 | _ -> assert(false) (* atom should have been skipped before *)

let read_name input =
  input_line input

let skip_header_lines input =
  let (_: string) = input_line input in
  let (_: string) = input_line input in
  ()

let read_atom_bonds_header input =
  let to_parse = input_line input in
  (* first two integers are 3-char fixed width each *)
  let num_atoms = int_of_string (S.strip (S.sub to_parse 0 3)) in
  let num_bonds = int_of_string (S.strip (S.sub to_parse 3 3)) in
  assert(S.ends_with to_parse "V2000");
  (num_atoms, num_bonds)

let read_bond_line input =
  let to_parse = input_line input in
  (* first two integers are 3-char fixed width each *)
  let src = int_of_string (S.strip (S.sub to_parse 0 3)) in
  let dst = int_of_string (S.strip (S.sub to_parse 3 3)) in
  (* indexes start at 1 in the SDF but atom indexes in the atoms array
   * start at 0 *)
  (src - 1, dst - 1)

let parse_atom_line input =
  let to_parse = input_line input in
  (* "^    5.0751   -3.8284   -4.0739 Br  0  0  0  0  0  0  0  0  0  0  0  0$" *)
  try
    Scanf.sscanf to_parse
      " %f %f %f %s@ " (* ignore all the rest of the line *)
      (fun x y z elt_symbol ->
         (anum_of_symbol elt_symbol, Vector3.make x y z))
  with exn ->
    let () = Log.fatal "Sdf_3D.parse_atom_line: cannot parse: '%s'"
        to_parse in
    raise exn

exception Four_dollars

let read_one_molecule input =
  let name = read_name input in
  (skip_header_lines input);
  let num_atoms, num_bonds = read_atom_bonds_header input in
  let elements = Array.make num_atoms 0 in
  let coords =
    Array.init num_atoms
      (fun i ->
         let (anum, xyz) = parse_atom_line input in
         elements.(i) <- anum;
         xyz
      ) in
  (* read all bonds *)
  let bonds = A.create num_atoms [] in
  for _i = 1 to num_bonds do
    let src, dst = read_bond_line input in
    bonds.(src) <- dst :: bonds.(src)
  done;
  (* put them back in the same order than what was read *)
  for i = 0 to num_atoms - 1 do
    bonds.(i) <- L.rev bonds.(i)
  done;
  (try
     (* look for end of this molecule's record *)
     while true do
       if input_line input = "$$$$" then
         raise Four_dollars
     done;
     assert(false)
   with Four_dollars -> ()
  );
  { name; elements; coords; bonds }

(* find all atoms within cutoff distance of the center atom *)
(* WARNING: this has O(n^2) complexity; we might index the atoms
 *   into some computational geometry data structure later on to accelerate
 * this if bottleneck *)
let within_cutoff cut mol i_atom =
  let cut2 = cut *. cut in
  let center = mol.coords.(i_atom) in
  A.fold_lefti (fun acc i coord ->
      if i = i_atom then
        (* also include the atom itself: since there must
         * be one channel with an atom at distance
         * 0.0 to distinguish different chemical elements being
         * atom centers *)
        let anum = mol.elements.(i) in
        (i, anum, coord) :: acc
      else
        begin
          let dist2 = V3.mag2 (V3.diff center coord) in
          if dist2 < cut2 then
            let anum = mol.elements.(i) in
            (i, anum, coord) :: acc
          else
            acc
        end
    ) [] mol.coords

let connected_atoms mol i_atom =
  let connected = mol.bonds.(i_atom) in
  L.map (fun dst ->
      let anum = mol.elements.(dst) in
      let coord = mol.coords.(dst) in
      (dst, anum, coord)
    ) connected

(* FBR: possible evolution 3: use a vanishing kernel to weight contributions
        instead of just a hard cutoff distance *)

let encode_first_layer dx cutoff mol =
  let nx = 1 + int_of_float (cutoff /. dx) in
  (* Log.debug "nx: %d" nx; *)
  let nb_atoms = A.length mol.elements in
  (* just encode the center atom's radial environment *)
  A.init nb_atoms (fun atom_i ->
      let res = A.make_matrix nx nb_channels 0.0 in
      let center = mol.coords.(atom_i) in
      let neighbors = within_cutoff cutoff mol atom_i in
      L.iter (fun (_atom_j, anum, coord) ->
          if anum < 0 then
            () (* unsupported elt. already reported before *)
            else
              let chan = channel_of_anum anum in
              (* Log.debug "chan: %d" chan; *)
              let dist = V3.dist center coord in
              (* Log.debug "dist: %g" dist; *)
              let bin_before = int_of_float (dist /. dx) in
              (* Log.debug "x_l: %d" bin_before; *)
              let bin_after = bin_before + 1 in
              (* Log.debug "x_r: %d" bin_after; *)
              let before = dx *. (float bin_before) in
              (* Log.debug "before: %g" before; *)
              let after = before +. dx in
              (* Log.debug "after: %g" after; *)
              (* linear binning *)
              let w_l = 1.0 -. (dist -. before) in
              let w_r = 1.0 -. (after -. dist) in
              res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
              res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
        ) neighbors;
      res
    )

(* evolution 2: as a second layer, add envs. from all neighbor
   atoms (more crowded than evol 1) *)
let encode_two_layers dx cutoff mol =
  let nx = 1 + int_of_float (cutoff /. dx) in
  let nb_atoms = A.length mol.elements in
  A.init nb_atoms (fun atom_i ->
      let res = A.make_matrix nx (2 * nb_channels) 0.0 in
      let center = mol.coords.(atom_i) in
      let center' = (atom_i, mol.elements.(atom_i), center) in
      let neighbors = within_cutoff cutoff mol atom_i in
      (* encode the center atom's radial environment *)
      L.iter (fun (_atom_j, anum, coord) ->
          if anum < 0 then
            () (* unsupported elt. already reported before *)
            else
              let chan = channel_of_anum anum in
              let dist = V3.dist center coord in
              let bin_before = int_of_float (dist /. dx) in
              let bin_after = bin_before + 1 in
              let before = dx *. (float bin_before) in
              let after = before +. dx in
              (* linear binning *)
              let w_l = 1.0 -. (dist -. before) in
              let w_r = 1.0 -. (after -. dist) in
              res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
              res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
        ) neighbors;
      (* also encode environments of all its Cartesian neighbors;
         but only those within the sphere centered on the center atom *)
      let neighbors' = center' :: neighbors in
      L.iter (fun (atom_j, _anum_j, coord_j) ->
          let neighbors'' =
            L.filter (fun (atom_k, _anum_k, coord_k) ->
                (atom_j <> atom_k) && (* not same atom *)
                (V3.dist coord_j coord_k < cutoff) (* within cutoff distance *)
              ) neighbors' in
          (* encode this atom's radial environment *)
          L.iter (fun (_atom_l, anum_l, coord_l) ->
              if anum_l < 0 then
                () (* unsupported elt. already reported before *)
              else
                let chan = nb_channels + (channel_of_anum anum_l) in
                let dist = V3.dist coord_j coord_l in
                let bin_before = int_of_float (dist /. dx) in
                let bin_after = bin_before + 1 in
                let before = dx *. (float bin_before) in
                let after = before +. dx in
                (* linear binning *)
                let w_l = 1.0 -. (dist -. before) in
                let w_r = 1.0 -. (after -. dist) in
                res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
                res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
            ) neighbors''
        ) neighbors;
      res
    )

(* evolution 1: as a second layer, add envs. from all connected atoms *)
let encode_two_layers' dx cutoff mol =
  let nx = 1 + int_of_float (cutoff /. dx) in
  let nb_atoms = A.length mol.elements in
  A.init nb_atoms (fun atom_i ->
      let res = A.make_matrix nx (2 * nb_channels) 0.0 in
      let center = mol.coords.(atom_i) in
      let center' = (atom_i, mol.elements.(atom_i), center) in
      let neighbors = within_cutoff cutoff mol atom_i in
      (* encode the center atom's radial environment *)
      L.iter (fun (_atom_j, anum, coord) ->
          if anum < 0 then
            () (* unsupported elt. already reported before *)
            else
              let chan = channel_of_anum anum in
              let dist = V3.dist center coord in
              let bin_before = int_of_float (dist /. dx) in
              let bin_after = bin_before + 1 in
              let before = dx *. (float bin_before) in
              let after = before +. dx in
              (* linear binning *)
              let w_l = 1.0 -. (dist -. before) in
              let w_r = 1.0 -. (after -. dist) in
              res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
              res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
        ) neighbors;
      (* also encode environments of all its connected neighbors
         within the sphere centered on the center atom *)
      let neighbors' = connected_atoms mol atom_i in
      L.iter (fun (_atom_j, _anum_j, coord_j) ->
          (* encode this atom's radial environment *)
          L.iter (fun (_atom_l, anum_l, coord_l) ->
              if anum_l < 0 then
                () (* unsupported elt. already reported before *)
              else (* within cutoff distance *)
              if V3.dist coord_j coord_l < cutoff then
                let chan = nb_channels + (channel_of_anum anum_l) in
                let dist = V3.dist coord_j coord_l in
                let bin_before = int_of_float (dist /. dx) in
                let bin_after = bin_before + 1 in
                let before = dx *. (float bin_before) in
                let after = before +. dx in
                (* linear binning *)
                let w_l = 1.0 -. (dist -. before) in
                let w_r = 1.0 -. (after -. dist) in
                res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
                res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
            ) (center' :: neighbors')
        ) neighbors';
      res
    )

(* [nb_layers]: how far we should consider from the atom center
   (distance in bonds on the molecular graph).
   [cutoff]: how far from the center atom in Angstrom we should consider
   (in Cartesian space)
   [dx]: axis discretization step
   [mol]: molecule to encode *)
let encode_atoms (nb_layers: int) (cutoff: float) (dx: float) (mol: atoms_3D)
  : encoded_atom array =
  match nb_layers with
  | 1 -> encode_first_layer dx cutoff mol
  | 2 -> encode_two_layers dx cutoff mol
  | 3 -> encode_two_layers' dx cutoff mol
  | _ -> failwith "not implemented yet"
