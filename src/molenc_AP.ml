(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * New (2023) counted atom pairs encoder w/ optional molecular standardization. *)

open Printf

module CLI = Minicli.CLI
module LO = Line_oriented
module Log = Dolog.Log

type mode = Input_dict of string
          | Output_dict of string

let verbose = ref false

let molenc_std_py =
  let exit_code, path = BatUnix.run_and_read "which molenc_std.py" in
  let res = BatString.strip path in
  Log.info "found molenc_std.py at %s" res;
  if exit_code <> BatUnix.WEXITED 0 then
    failwith "molenc_AP: no molenc_std.py in PATH"
  else
    res

let standardize_molecules in_fn out_fn =
  let cmd = sprintf "%s -i %s -o %s" molenc_std_py in_fn out_fn in
  (if !verbose then Log.debug "running: %s" cmd);
  let exit_code = Unix.system cmd in
  if exit_code <> BatUnix.WEXITED 0 then
    Log.warn "Molenc_AP.standardize_molecules: error while running: %s" cmd

(* read [chunk_size] molecules and store them in a temp_file *)
let read_some chunk_size input =
  let count = ref 0 in
  let tmp_smi_fn = Filename.temp_file "" (* no_prefix *) ".smi" in
  LO.with_out_file tmp_smi_fn (fun output ->
      try
        for _i = 1 to chunk_size do
          let line = input_line input in
          fprintf output "%s\n" line;
          incr count
        done
      with End_of_file -> ()
    );
  if !count = 0 then
    (Sys.remove tmp_smi_fn;
     raise Parany.End_of_input)
  else
    tmp_smi_fn

let standardize_some tmp_smi_fn =
  let tmp_std_smi_fn = Filename.temp_file "" (* no_prefix *) "_std.smi" in
  standardize_molecules tmp_smi_fn tmp_std_smi_fn;
  Sys.remove tmp_smi_fn;
  tmp_std_smi_fn

let do_not_standardize_some tmp_smi_fn =
  tmp_smi_fn

let catenate_some dst_fn tmp_src_fn =
  let cmd = sprintf "cat %s >> %s" tmp_src_fn dst_fn in
  (if !verbose then Log.debug "running: %s" cmd);
  let exit_code = Unix.system cmd in
  (if exit_code <> BatUnix.WEXITED 0 then
     Log.warn "Molenc_AP.catenate_some: error while running: %s" cmd
  );
  Sys.remove tmp_src_fn

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  (if argc = 1 then
     (eprintf "usage:\n  \
               %s -i in.smi -o out.AP\n  \
               -i <input.smi>: input molecules\n  \
               -o <output.AP>: unfolded counted atom pairs output\n  \
               -d <dico.dix>: use existing feature dictionary\n  \
               [--no-std]: do not standardize molecules\n  \
               [-m <int>]: maximum atom pairs distance (in bonds; default=OFF)\n  \
               [-np <int>]: parallelize on NCORES (default=1)\n  \
               [-cs <int>]: chunk size (default=50)\n  \
               [-v]: verbose/debug mode\n"
        Sys.argv.(0);
      exit 1)
  );
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-cs"] args 50 in
  let no_std = CLI.get_set_bool ["--no-std"] args in
  let _max_dist_opt = CLI.get_int_opt ["-m"] args in
  let _dict_mode = match CLI.get_string_opt ["-d"] args with
    | None -> Output_dict (input_fn ^ ".dix")
    | Some fn -> Input_dict fn in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let standardization =
    if no_std then
      do_not_standardize_some
    else
      standardize_some in
  LO.with_in_file input_fn (fun input ->
      Parany.run ~preserve:true nprocs
        ~demux:(fun () -> read_some csize input)
        ~work:standardization
        ~mux:(catenate_some output_fn)
    )

let () = main ()
