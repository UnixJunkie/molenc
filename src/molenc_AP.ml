(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * New (Nov 2023) counted atom pairs encoder w/ optional molecular standardization. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IMap = BatMap.Int
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

type mode = Input_dict of string
          | Output_dict of string

type atom = int array

module Atom_pair = struct

  type t = { left: atom ;
             dist: int ; (* topological distance in bonds *)
             right: atom }

  let compare = compare

  let create t1 dist t2 =
    (* canonicalization here *)
    if t1 <= t2 then
      { left = t1; dist; right = t2 }
    else
      { left = t2; dist; right = t1 }

end

module APM = Map.Make (Atom_pair)

type unfold_count_fp = { name: string;
                         feat_counts: int APM.t }

let fp_string_output mode out dict fp =
  fprintf out "%s,0.0,[" fp.name;
  let feat_counts = match mode with
    | Output_dict _ -> (* writable dict *)
      (* feature index 0 is a reserved value for later users of the same dict:
         unkown feature *)
      let feat_i = ref (1 + Ht.length dict) in
      APM.fold (fun feat count acc ->
          let feat' =
            try Ht.find dict feat
            with Not_found ->
              (Ht.add dict feat !feat_i;
               incr feat_i;
               !feat_i - 1) in
          IMap.add feat' count acc
        ) fp.feat_counts IMap.empty
    | Input_dict _ -> (* read-only dict *)
      APM.fold (fun feat count acc ->
          try IMap.add (Ht.find dict feat) count acc
          with Not_found ->
            let () = Log.warn "unknown feat" in
            let prev_count =
              try IMap.find 0 acc
              with Not_found -> 0 in
            IMap.add 0 (1 + prev_count) acc
        ) fp.feat_counts IMap.empty in
  let started = ref false in
  IMap.iter (fun feat count ->
      if !started then
        fprintf out ";%d:%d" feat count
      else
        (started := true;
         fprintf out "%d:%d" feat count)
    ) feat_counts;
  fprintf out "]\n"

(* unfolded counted atom pairs fingerprint encoding *)
let encode_smiles_line line =
  let smi, name = BatString.split ~by:"\t" line in
  let mol = Rdkit.__init__ ~smi () in
  let n = Rdkit.get_num_atoms mol () in
  let atom_types = A.init n (fun i -> Rdkit.type_EltFCaroNeighbs mol ~i ()) in
  let fp = ref APM.empty in
  (* autocorrelation *)
  for i = 0 to n - 2 do
    for j = i + 1 to n - 1 do
      let feat =
        Atom_pair.create
          atom_types.(i) (Rdkit.get_distance mol ~i ~j ()) atom_types.(j) in
      fp :=
        try APM.add feat (1 + APM.find feat !fp) !fp
        with Not_found -> APM.add feat 1 !fp
    done
  done;
  { name; feat_counts = !fp }

let verbose = ref false

(* find where is installed the molenc_std.py script *)
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
    (assert(0 = Sys.command (sprintf "rm -rf %s" tmp_dir));
     raise Parany.End_of_input)
  else
    tmp_dir

let standardize_some do_not tmp_dir =
  let smi_fn     = sprintf "%s/in.smi"     tmp_dir in
  let std_smi_fn = sprintf "%s/in_std.smi" tmp_dir in
  (if do_not then
     assert(0 = Sys.command (sprintf "mv %s %s" smi_fn std_smi_fn))
   else
     standardize_molecules smi_fn std_smi_fn
  );
  tmp_dir

(* molecular encoding, but feature dictionary does not exist yet *)
let encode_some tmp_dir =
  let smi_fn = sprintf "%s/in_std.smi" tmp_dir in
  let apb_fn = sprintf "%s/in_std.apb" tmp_dir in
  LO.save apb_fn (LO.map smi_fn encode_smiles_line)

let catenate_some dst_fn tmp_dir =
  (* FBR: CHANGE HERE ONCE ENCODER READY *)
  let cmd = sprintf "cat %s/in_std.smi >> %s" tmp_dir dst_fn in
  (if !verbose then Log.debug "running: %s" cmd);
  let exit_code = Unix.system cmd in
  (if exit_code <> BatUnix.WEXITED 0 then
     Log.warn "Molenc_AP.catenate_some: error while running: %s" cmd
  );
  assert(0 = Sys.command (sprintf "rm -rf %s" tmp_dir))

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
               [-f]: overwrite output file, if any\n  \
               [--no-std]: do not standardize molecules\n  \
               [-m <int>]: maximum atom pairs distance (in bonds; default=OFF)\n  \
               [-np <int>]: parallelize on NCORES (default=1)\n  \
               [-c <int>]: chunk size (default=50)\n  \
               [-v]: verbose/debug mode\n"
        Sys.argv.(0);
      exit 1)
  );
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 50 in
  let no_std = CLI.get_set_bool ["--no-std"] args in
  let force = CLI.get_set_bool ["-f"] args in
  let _max_dist_opt = CLI.get_int_opt ["-m"] args in
  let _dict_mode = match CLI.get_string_opt ["-d"] args with
    | None -> Output_dict (input_fn ^ ".dix")
    | Some fn -> Input_dict fn in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (if Sys.file_exists output_fn then
     let () = Log.warn "Molenc_AP.main: output file exists: %s" output_fn in
     if not force then
       exit 1
     else
       Sys.remove output_fn
  );
  LO.with_in_file input_fn (fun input ->
      Parany.run ~preserve:true nprocs
        ~demux:(fun () -> read_some csize input)
        ~work:(standardize_some no_std)
        ~mux:(catenate_some output_fn)
    )

let () = main ()
