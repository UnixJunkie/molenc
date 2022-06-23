(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Read molecule name, element symbols and their 3D coordinates from a .sdf file
   holding 3D conformers. *)

open Printf

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

type encoded_atom =
  { radial: float array array; (* shape: (nb_dx, nb_chans) *)
    angular: float array array (* shape: (nb_da, nb_chans +
                                  nb_chans*(nb_chans-1)/2 *) }

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
  | x ->
    let () = Log.warn "anum: %d" x in
    10

let nb_angular_channels =
  nb_channels + (nb_channels * (nb_channels - 1) / 2)

let angular_channel_of_anums a1 a2 =
  let (x, y) =
    if a1 <= a2 then (a1, a2)
    else (a2, a1) in
  match (x, y) with
  | (1 , 1 ) -> 0
  | (6 , 6 ) -> 1
  | (7 , 7 ) -> 2
  | (8 , 8 ) -> 3
  | (9 , 9 ) -> 4
  | (15, 15) -> 5
  | (16, 16) -> 6
  | (17, 17) -> 7
  | (35, 35) -> 8
  | (53, 53) -> 9
  | (1 , 6 ) -> 10
  | (1 , 7 ) -> 11
  | (1 , 8 ) -> 12
  | (1 , 9 ) -> 13
  | (1 , 15) -> 14
  | (1 , 16) -> 15
  | (1 , 17) -> 16
  | (1 , 35) -> 17
  | (1 , 53) -> 18
  | (6 , 7 ) -> 19
  | (6 , 8 ) -> 20
  | (6 , 9 ) -> 21
  | (6 , 15) -> 22
  | (6 , 16) -> 23
  | (6 , 17) -> 24
  | (6 , 35) -> 25
  | (6 , 53) -> 26
  | (7 , 8 ) -> 27
  | (7 , 9 ) -> 28
  | (7 , 15) -> 29
  | (7 , 16) -> 30
  | (7 , 17) -> 31
  | (7 , 35) -> 32
  | (7 , 53) -> 33
  | (8 , 9 ) -> 34
  | (8 , 15) -> 35
  | (8 , 16) -> 36
  | (8 , 17) -> 37
  | (8 , 35) -> 38
  | (8 , 53) -> 39
  | (9 , 15) -> 40
  | (9 , 16) -> 41
  | (9 , 17) -> 42
  | (9 , 35) -> 43
  | (9 , 53) -> 44
  | (15, 16) -> 45
  | (15, 17) -> 46
  | (15, 35) -> 47
  | (15, 53) -> 48
  | (16, 17) -> 49
  | (16, 35) -> 50
  | (16, 53) -> 51
  | (17, 35) -> 52
  | (17, 53) -> 53
  | (35, 53) -> 54
  | (_ , _ ) -> assert(false)

let symbol_of_channel = [|"C";"H";"N";"O";"P";"S";"F";"Cl";"Br";"I"|]

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

(* like [within_cutoff]; but the i_atom is not included in the list *)
let within_cutoff' cut mol i_atom =
  let cut2 = cut *. cut in
  let center = mol.coords.(i_atom) in
  A.fold_lefti (fun acc i coord ->
      if i = i_atom then
        acc
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

let pow3 x =
  x *. x *. x

let pow2 x =
  x *. x

let triweight_K x =
  if x < 1.0 then
    (* 1.09375 = 35/32 *)
    1.09375 *. pow3 (1.0 -. (pow2 x))
  else
    assert(false)

(* evaluate normalized kernel *)
let eval_K bwidth x =
  if x >= bwidth then
    0.0
  else
    let scale = 1.0 /. bwidth in
    let x' = x /. bwidth in
    scale *. (triweight_K x')

(* FBR: output encoding to a data file for gnuplot *)

let pi = 4.0 *. atan 1.0

(* compute all angles involving given center atom
   and all its spatial neighbors *)
let all_angles (center: V3.t) (neighbors: (int * int * V3.t) list)
  : (int * float) list =
  let angles (_, _, x_xyz) others =
    let xc = V3.diff x_xyz center in
    L.rev_map (fun (_, z_anum, z_xyz) ->
        let zc = V3.diff z_xyz center in
        let angle = V3.angle xc zc in
        assert(0.0 <= angle && angle <= pi);
        (z_anum, angle)
      ) others in
  let rec loop acc = function
    | [] -> acc
    | x :: xs ->
      let ys = angles x xs in
      let acc' = L.rev_append ys acc in
      loop acc' xs in
  loop [] neighbors

let encode_first_layer dx cutoff da mol =
  let nx = 1 + int_of_float (cutoff /. dx) in
  let na = 1 + int_of_float (pi /. da) in
  (* Log.debug "nx: %d" nx; *)
  let nb_atoms = A.length mol.elements in
  (* encode radial environment around center atom *)
  let radial_envs =
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
      ) in
  (* encode all angles involving the center atom *)
  let angular_envs =
    A.init nb_atoms (fun atom_i ->
        let res = A.make_matrix na nb_angular_channels 0.0 in
        let center = mol.coords.(atom_i) in
        let center_anum = mol.elements.(atom_i) in
        let neighbors = within_cutoff' cutoff mol atom_i in
        let angles = all_angles center neighbors in
        L.iter (fun (anum, angle) ->
            if anum < 0 then
              () (* unsupported elt. already reported before *)
            else
              let chan = angular_channel_of_anums center_anum anum in
              (* Log.debug "chan: %d" chan; *)
              let bin_before = int_of_float (angle /. da) in
              (* Log.debug "x_l: %d" bin_before; *)
              let bin_after = bin_before + 1 in
              (* Log.debug "x_r: %d" bin_after; *)
              let before = da *. (float bin_before) in
              (* Log.debug "before: %g" before; *)
              let after = before +. da in
              (* Log.debug "after: %g" after; *)
              (* linear binning *)
              let w_l = 1.0 -. (angle -. before) in
              let w_r = 1.0 -. (after -. angle) in
              res.(bin_before).(chan) <- res.(bin_before).(chan) +. w_l;
              res.(bin_after).(chan) <- res.(bin_after).(chan) +. w_r
          ) angles;
        res
      ) in
  A.map2 (fun radial angular ->
      { radial; angular }
    ) radial_envs angular_envs

(* [nb_layers]: how far we should consider from the atom center
   (distance in bonds on the molecular graph).
   [cutoff]: how far from the center atom in Angstrom we should consider
   (in Cartesian space)
   [dx]: axis discretization step
   [mol]: molecule to encode *)
let encode_atoms
    (nb_layers: int) (cutoff: float) (dx: float) (da: float) (mol: atoms_3D)
  : encoded_atom array =
  match nb_layers with
  | 1 -> encode_first_layer dx cutoff da mol
  | _ -> failwith (sprintf "unsupported l: %d" nb_layers)
