(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Dump a .sdf file in txt format. *)

open Printf

module A = BatArray
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
          printf "%s\n" name;
          A.iter2 (fun xyz anum ->
              let (x, y, z) = V3.to_triplet xyz in
              let elt = Sdf_3D.symbol_of_anum anum in
              printf "%f %f %f %s\n" x y z elt
            ) coords elts
        done
      with End_of_file -> ()
    )

let () = main ()
