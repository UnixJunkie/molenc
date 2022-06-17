(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Read molecule name, element symbols and their 3D coordinates from a .sdf file
   holding 3D conformers. *)

module Log = Dolog.Log
module S = BatString

type t = { name: string;
           elements: int array;
           coords: Vector3.t array }

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
  | _ -> -1 (* unsupported elt. *)

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
  (* first integer is 3-char fixed width *)
  let num_atoms = int_of_string (S.strip (S.left to_parse 3)) in
  let num_bonds = int_of_string (S.strip (S.left (S.lchop ~n:3 to_parse) 3)) in
  assert(S.ends_with to_parse "V2000");
  (num_atoms, num_bonds)

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
  (* skip all bonds *)
  for _i = 1 to num_bonds do
    ignore(input_line input)
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
  { name; elements; coords }
