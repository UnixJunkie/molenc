(* Copyright (C) 2022, Francois Berenger

   Tsuda laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Dump a 3D .sdf file's atoms using 3DAE. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Sdf_3D = Molenc.Sdf_3D

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -i <filename.sdf>: input file\n  \
              -o <filename.csv>: output file\n  \
              [--charges <charges.csv>]: prepend those charges\n  \
              [--separate]: separate channels in output files\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-l <int>]: number of layers in [1,2]\n  \
              [-c <float>]: cutoff distance (default=3.5A)\n  \
              [-dx <float>]: discretization step (default=0.1A)\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let maybe_charges_fn = CLI.get_string_opt ["--charges"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let nb_layers = CLI.get_int_def ["-l"] args 1 in
  let cutoff = CLI.get_float_def ["-c"] args 3.5 in
  let dx = CLI.get_float_def ["-dx"] args 0.1 in
  let separate_channels = CLI.get_set_bool ["--separate"] args in
  let max_feat = ref (-1) in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let charges = match maybe_charges_fn with
    | None -> [||]
    | Some charges_fn ->
      let lines = LO.lines_of_file charges_fn in
      let n = L.length lines in
      Log.info "%s: %d charges" charges_fn n;
      let res = A.create n 0.0 in
      L.iteri (fun i charge_str ->
          res.(i) <- float_of_string charge_str
        ) lines;
      res in
  let per_channel_output_files =
    if separate_channels then
      A.init Sdf_3D.nb_channels (fun chan ->
          let symbol = Sdf_3D.symbol_of_channel.(chan) in
          let chan_out_fn = sprintf "%s.%s" output_fn symbol in
          open_out chan_out_fn
        )
    else [||] in
  let dev_null_chan = open_out "/dev/null" in
  let prepend_charges = A.length charges > 0 in
  let atoms_count = ref 0 in
  LO.with_infile_outfile input_fn output_fn (fun input output' ->
      try
        while true do
          let mol = Sdf_3D.read_one_molecule input in
          let atoms_3dae = Sdf_3D.encode_atoms nb_layers cutoff dx mol in
          A.iteri (fun i_atom encoded_atom ->
              let radial = Sdf_3D.(encoded_atom.radial) in
              let output =
                if separate_channels then
                  begin
                    let anum = mol.elements.(i_atom) in
                    let chan = Sdf_3D.channel_of_anum anum in
                    if chan < Sdf_3D.nb_channels then
                      per_channel_output_files.(chan)
                    else
                      dev_null_chan
                  end
                else
                  output' in
              (if prepend_charges then
                 let () = fprintf output "%f" charges.(!atoms_count) in
                 incr atoms_count
               else
                 fprintf output "%d" i_atom
              );
              let nb_dx = A.length radial in
              let nb_chans = A.length radial.(0) in
              for i_chan = 0 to nb_chans - 1 do
                for i_dx = 0 to nb_dx - 1 do
                  let feat = radial.(i_dx).(i_chan) in
                  if feat > 0.0 then
                    (* the feature vector should be very sparse;
                         liblinear wants feature indexes to start at 1 *)
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
  Log.info "num_atoms: %d max_feat: %d" !atoms_count !max_feat;
  (* close separate channels, if any *)
  A.iter close_out per_channel_output_files;
  close_out dev_null_chan

let () = main ()
