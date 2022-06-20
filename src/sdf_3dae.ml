(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Dump a 3D .sdf file's atoms using 3DAE. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module LO = Line_oriented
module Log = Dolog.Log
module Sdf_3D = Molenc.Sdf_3D
module V3 = Vector3

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename.sdf>]: input file\n  \
              -o <filename.csv>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-l <int>]: number of layers in [1,2]\n  \
              [-c <float>]: cutoff distance (default=3.5A)\n  \
              [-dx <float>]: discretization step (default=0.1A)\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let output_fn = CLI.get_string ["-o"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let nb_layers = CLI.get_int_def ["-l"] args 1 in
  let cutoff = CLI.get_float_def ["-c"] args 3.5 in
  let dx = CLI.get_float_def ["-dx"] args 0.1 in
  let max_feat = ref (-1) in
  CLI.finalize (); (* ------------------------------------------------------ *)
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let mol = Sdf_3D.read_one_molecule input in
          let atoms_3dae = Sdf_3D.encode_atoms nb_layers cutoff dx mol in
          A.iteri (fun i_atom encoded_atom ->
              fprintf output "%d" i_atom;
              let nb_dx = A.length encoded_atom in
              let nb_chans = A.length encoded_atom.(0) in
              for i_chan = 0 to nb_chans - 1 do
                for i_dx = 0 to nb_dx - 1 do
                  let feat = encoded_atom.(i_dx).(i_chan) in
                  if feat > 0.0 then
                    (* the feature vector should be very sparse;
                       lilinear wants feature indexes to start at 1 *)
                    let feat_idx = 1 + i_dx + (i_chan * nb_dx) in
                    (if feat_idx > !max_feat then
                       max_feat := feat_idx
                    );
                    fprintf output " %d:%g" feat_idx feat
                done
              done;
              fprintf output "\n"
            ) atoms_3dae
        done
      with End_of_file -> ()
    );
  Log.info "max feature index: %d" !max_feat

let () = main ()
