(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Dump a .sdf file in txt format. *)

open Printf

module A = BatArray
module L = BatList
module LO = Line_oriented
module Sdf_3D = Molenc.Sdf_3D
module V3 = Vector3

let main () =
  let input_fn = Sys.argv.(1) in
  LO.with_in_file input_fn (fun input ->
      try
        while true do
          let mol = Sdf_3D.read_one_molecule input in
          let name = Sdf_3D.(mol.name) in
          let elts = Sdf_3D.(mol.elements) in
          let coords = Sdf_3D.(mol.coords) in
          let bonds = Sdf_3D.(mol.bonds) in
          printf "%s\n" name;
          A.iter2 (fun xyz anum ->
              let (x, y, z) = V3.to_triplet xyz in
              let elt = Sdf_3D.symbol_of_anum anum in
              printf "%10.4f%10.4f%10.4f %s\n" x y z elt
            ) coords elts;
          A.iteri (fun src_a connected_atoms ->
              L.iter (fun dst_a ->
                  printf "%3d%3d\n" (1 + src_a) (1 + dst_a)
                ) connected_atoms
            ) bonds
        done
      with End_of_file -> ()
    )

let () = main ()
