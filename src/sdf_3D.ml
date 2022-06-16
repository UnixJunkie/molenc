(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Read molecule name, element symbols and their 3D coordinates from a .sdf file
   holding 3D conformers. *)

type t = { name: string;
           coords: Vector3.t array }

(* FBR: add a driver program for tests: molenc_sdf_read *)

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
  with exn -> raise exn
