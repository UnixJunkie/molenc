(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * RFP encoder *)

(* FBR: carefully test on some molecules
   - use an encoding dictionary; reserve feature 0 for unknown features;
     maybe integrate into molenc_AP; w/ CLI option like --UCECFP
 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IMap = BatMap.Int
module SMap = BatMap.String
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString
module Utls = Molenc.Utls

(* the Rdkit module relies on Pyml *)
let () = Py.initialize ~version:3 ()

type mode = Input
          | Output

type fp = { name: string;
            feat_counts: int SMap.t }

(* chemical formula at [radius] bonds away from [center_atom_i] *)
let get_atom_env distances radius indexes elements =
  A.fold (fun acc a_i ->
      if distances.(a_i) = radius then
        let elt = A.unsafe_get elements a_i in
        let prev_count = SMap.find_default 0 elt acc in
        SMap.add elt (prev_count + 1) acc
      else (* not part of the atom_env *)
        acc
    ) SMap.empty indexes

(* chemical formula for atom env. layer *)
let string_of_elt2count smap =
  let buff = Buffer.create 128 in
  SMap.iter (fun elt count ->
      if count > 1 then
        bprintf buff "%s%d" elt count
      else (* count=1 *)
        Buffer.add_string buff elt
    ) smap;
  Buffer.contents buff

(* algorithm to encode given atom:
   - get max radius away from this atom (R)
   - trim to max_radius the user is asking (M)
   - get all layers of atom environments for this atom as a list
   - record features concerning this atom by varying the radius
     from 0 to min(R,M)

   to encode the whole molecule, we need to do this for each atom

   then we count the number of times each atom env. was seen
 *)
let environments_for_atom distances fp_radius elements =
  let max_radius = min fp_radius (A.max distances) in
  (* r=0 must also be taken into account *)
  let res = A.init (1 + max_radius) (fun _rad -> SMap.empty) in
  A.iteri (fun i dist ->
      if dist <= max_radius then
        let elt = A.unsafe_get elements i in
        let smap = A.unsafe_get res dist in
        let count = SMap.find_default 0 elt smap in
        A.unsafe_set res dist (SMap.add elt (count + 1) smap)
    ) distances;
  let formulas = A.to_list (A.map string_of_elt2count res) in
  (* return all strings; from r=0 to max_radius *)
  let l = ref [] in
  for i = 1 + max_radius downto 1 do
    l := (String.concat "," (L.take i formulas)) :: !l
  done;
  !l

let environments_for_molecule fp_radius mol_noH =
  let mol = Rdkit.add_hydrogens mol_noH () in
  let num_atoms = Rdkit.get_num_atoms mol () in
  let elements = Rdkit.get_elements mol () in
  let atom_envs =
    (* gather all environments *)
    A.init num_atoms (fun i ->
        let dists = Rdkit.get_distances mol ~i () in
        environments_for_atom dists fp_radius elements
      ) in
  (* count them *)  
  A.fold_left (fun acc0 l ->
      L.fold_left (fun acc1 formula ->
          let count = SMap.find_default 0 formula acc1 in
          SMap.add formula (count + 1) acc1
        ) acc0 l
    ) SMap.empty atom_envs

(* unfolded counted RFP encoding *)
let encode_smiles_line max_radius line =
  let smi, name = BatString.split ~by:"\t" line in
  let mol_noH = Rdkit.__init__ ~smi () in
  { name; feat_counts = environments_for_molecule max_radius mol_noH }


let fp_string_output mode out dict fp =
  fprintf out "%s,0.0,[" fp.name;
  let feat_counts = match mode with
    | Output -> (* writable dict *)
      (* feature index 0 is a reserved value for later users of the same dict:
         unkown feature *)
      let feat_i = ref (1 + Ht.length dict) in
      SMap.fold (fun feat count acc ->
          let feat' =
            try Ht.find dict feat
            with Not_found ->
              let idx = !feat_i in
              Ht.add dict feat !feat_i;
              incr feat_i;
              idx in
          IMap.add feat' count acc
        ) fp.feat_counts IMap.empty
    | Input -> (* read-only dict *)
      SMap.fold (fun feat count acc ->
          try IMap.add (Ht.find dict feat) count acc
          with Not_found ->
            let () = Log.warn "unknown feat: %s" feat in
            let prev_count = IMap.find_default 0 0 acc in
            IMap.add 0 (1 + prev_count) acc
        ) fp.feat_counts IMap.empty in
  (* FBR: warn about the number of unknown features in a molecule *)  
  let started = ref false in
  IMap.iter (fun feat count ->
      if !started then
        fprintf out ";%d:%d" feat count
      else
        (started := true;
         fprintf out "%d:%d" feat count)
    ) feat_counts;
  fprintf out "]\n"

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
               [-o <output.fp>]: UCECFP output\n  \
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
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let _csize = CLI.get_int_def ["-c"] args 200 in
  let _force = CLI.get_set_bool ["-f"] args in
  let max_radius = match CLI.get_int_opt ["-m"] args with
    | None -> max_int
    | Some x -> x in
  let verbose = CLI.get_set_bool ["-v"] args in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (if verbose then
     Log.(set_log_level DEBUG)
  );
  (*
  let already_out_fn = Sys.file_exists output_fn in
  (if already_out_fn then
     (Log.warn "Molenc_AP: output file exists: %s" output_fn;
      if not force then
        exit 1
     );
  );
  *)
  let reads = ref 0 in
  let writes = ref 0 in
  let dict = Ht.create 20_011 in
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let line = input_line input in
          incr reads;
          let fp = encode_smiles_line max_radius line in
          fp_string_output Output output dict fp;
          incr writes
        done
      with End_of_file -> ()
    );
  let dt = Unix.gettimeofday() -. start in
  Log.info "wrote %#d molecules (%.2f Hz)" !writes ((float !writes) /. dt)

let () = main ()
