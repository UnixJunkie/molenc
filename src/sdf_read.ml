(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Dump a .sdf file in txt format. *)

module LO = Line_oriented

let main () =
  let input_fn = Sys.argv.(1) in
  LO.with_in_file input_fn (fun input ->
      try
        while true do
          let mol = Sdf_3D.read_one_molecule input in
          ()
        done
      with End_of_file -> ()
    )

let () = main ()
