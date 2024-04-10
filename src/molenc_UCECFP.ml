(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * UCECFP* encoder *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module SMap = BatMap.String
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit
module Utls = Molenc.Utls

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

type mode = Input
          | Output

type unfolded_counted_fp = { name: string;
                             feat_counts: int SMap.t }

(* FBR: implement bound on max atom env. radius *)

(* UECFP* encoding *)
let encode_smiles_line _max_radius line =
  let smi, _name = BatString.split ~by:"\t" line in
  let mol = Rdkit.__init__ ~smi () in
  let _num_atoms = Rdkit.get_num_atoms mol () in
  (* compute diameter *)
  (* encode each atom using all diameters from 0 to N *)
  failwith "not implemented yet"

let verbose = ref false

(* read [chunk_size] molecules and store them in a temp_file *)
let read_some reads chunk_size input =
  let count = ref 0 in
  let tmp_dir = Filename.temp_dir "molenc_" "" (* no suffix *) in
  let tmp_smi_fn = sprintf "%s/in.smi" tmp_dir in
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
    let () = Log.info "read %#d molecules" !reads in
    assert(0 = Sys.command (sprintf "rm -rf %s" tmp_dir));
    raise Parany.End_of_input
  else
    (reads := !reads + !count;
     tmp_dir)

(* molecular encoding, but feature dictionary does not exist yet *)
let encode_some max_dist simple_types tmp_dir =
  let smi_fn = sprintf "%s/in_std.smi" tmp_dir in
  LO.map smi_fn (encode_smiles_line max_dist simple_types)

let main () =
  let start = Unix.gettimeofday () in
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  (if argc = 1 then
     (eprintf "usage:\n  \
               %s -i in.smi -o out.fp\n  \
               -i <input.smi>: input molecules\n  \
               -o <output.fp>: UCECFP output\n  \
               [-f]: overwrite output file, if any\n  \
               [-m <int>]: maximum atom env. radius (in bonds; default=OFF)\n  \
               [-np <int>]: parallelize on NCORES (default=1)\n  \
               [-c <int>]: chunk size (default=200)\n  \
               [-v]: verbose/debug mode\n"
        Sys.argv.(0);
      exit 1)
  );
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 200 in
  let force = CLI.get_set_bool ["-f"] args in
  let max_radius = match CLI.get_int_opt ["-m"] args with
    | None -> max_int
    | Some x -> x in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let already_out_fn = Sys.file_exists output_fn in
  (if already_out_fn then
     begin
       Log.warn "Molenc_AP: output file exists: %s" output_fn;
       if not force then
         exit 1
     end
  );
  let reads = ref 0 in
  let writes = ref 0 in
  LO.with_infile_outfile input_fn output_fn (fun input _output ->
      (* !!! KEEP ~csize:1 below !!! *)
      Parany.run ~csize:1 ~preserve:true nprocs
        ~demux:(fun () -> read_some reads csize input)
        ~work:(fun tmp_dir -> encode_some max_radius tmp_dir)
        ~mux:(fun _ -> failwith "not implemented yet")
    );
  let dt = Unix.gettimeofday() -. start in
  Log.info "wrote %#d molecules (%.2f Hz)" !writes ((float !writes) /. dt)

let () = main ()
