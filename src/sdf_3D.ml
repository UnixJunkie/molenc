(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Read molecule name, element symbols and their 3D coordinates from a .sdf file
   holding 3D conformers. *)

module Log = Dolog.Log

type t = { name: string;
           elements: int array;
           coords: Vector3.t array }

(* FBR: add a driver program for tests: molenc_sdf_read *)

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

let read_name input =
  input_line input

let skip_header_lines input =
  let (_: string) = input_line input in
  let (_: string) = input_line input in
  ()

let read_atom_bonds_header input =
  let to_parse = input_line input in
  try Scanf.sscanf to_parse
        (* eg: "^ 29 30  0  0  0  0  0  0  0  0999 V2000$" *)
        "%d %d %d %d %d %d %d %d %d %d %s"
        (fun num_atoms num_bonds _ _ _ _ _ _ _ _ version ->
           assert(version = "V2000");
           (num_atoms, num_bonds)
        )
  with exn ->
    let () = Log.fatal "Sdf_3D.read_atom_bonds_header: cannot parse: %s"
        to_parse in
    raise exn

let parse_atom_line input =
  let to_parse = input_line input in
  (* "^    5.0751   -3.8284   -4.0739 Br  0  0  0  0  0  0  0  0  0  0  0  0$" *)
  try
    Scanf.sscanf to_parse
      "%f %f %f %s@ %d %d %d %d %d %d %d %d %d %d %d %d"
      (fun x y z elt_symbol _ _ _ _ _ _ _ _ _ _ _ _ ->
         (anum_of_symbol elt_symbol, Vector3.make x y z))
  with exn ->
    let () = Log.fatal "Sdf_3D.parse_atom_line: cannot parse: %s"
        to_parse in
    raise exn
