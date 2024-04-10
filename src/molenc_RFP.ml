(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * RFP encoder *)

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

(* FBR: atom env module? *)

(* chemical formula at [radius] away from [center_atom_i] *)
let get_atom_env center_atom_i radius mol indexes elements =
  A.fold (fun acc a_i ->
      (* FBR: export this distance matrix once and for all from Python *)
      if Rdkit.get_distance mol ~i:center_atom_i ~j:a_i () = radius then
        let elt = A.unsafe_get elements a_i in
        let prev_count = SMap.find_default 0 elt acc in
        SMap.add elt (prev_count + 1) acc
      else (* not part of the atom_env *)
        acc
    ) SMap.empty indexes

(* UECFP* encoding *)
let encode_smiles_line max_radius line =
  let smi, _name = BatString.split ~by:"\t" line in
  let mol = Rdkit.__init__ ~smi () in
  let num_atoms = Rdkit.get_num_atoms mol () in
  let diameter = Rdkit.get_diameter mol () in
  let elements = Rdkit.get_elements mol () in
  (* encode each atom using all diameters from 0 to max_radius *)
  let indexes = A.init num_atoms (fun i -> i) in
  let radii = A.init (min max_radius diameter) (fun i -> i) in
  let res = Buffer.create 1024 in
  SMap.iter (Printf.bprintf res "%s:%d")
    (A.fold (fun acc0 a_i ->
         let buff = Buffer.create 128 in
         A.fold (fun acc1 radius ->
             let atom_env = get_atom_env a_i radius mol indexes elements in
             (* get a chemical formula *)
             SMap.iter (Printf.bprintf buff "%s:%d") atom_env;
             let curr_env = Buffer.contents buff in
             let count = SMap.find_default 0 curr_env acc1 in
             SMap.add curr_env (count + 1) acc1
           ) acc0 radii
       ) SMap.empty indexes
    );
  Buffer.contents res

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
     (Log.warn "Molenc_AP: output file exists: %s" output_fn;
      if not force then
        exit 1
     );
  );
  let reads = ref 0 in
  let writes = ref 0 in
  let dt = Unix.gettimeofday() -. start in
  Log.info "wrote %#d molecules (%.2f Hz)" !writes ((float !writes) /. dt)

let () = main ()
