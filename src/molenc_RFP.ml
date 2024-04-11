(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * RFP encoder *)

(* FBR: in verbose mode: output atom environments by incr. (grouped by) radius *)

(* FBR: carefully test on some molecules *)

(* FBR: how to encode an atom environment as an unambiguous integer (and be reversible)?
        - I could assign a distinct prime number to each chemical element
        - then the sum of powers (the power is the number of times the element was seen)
          of those primes should be decodable
        - ask the CCL?
        - OR, just encode the corresp. string as a number?
        - Huffman compress the string, then corresponding integer?
        - create the right alphabet (all required chars plus 0-9 digits); then
          encode as a number in that base
 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module SMap = BatMap.String
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString
module Utls = Molenc.Utls

(* because the Rdkit module uses Pyml *)
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

let fp_to_string fp =
  let buff = Buffer.create 1024 in
  Printf.bprintf buff "%s\t" fp.name;
  let started = ref false in
  SMap.iter (fun env count ->
      (* separate environments *)
      if !started then
        Buffer.add_char buff ';'
      else
        started := true
      ;
      Printf.bprintf buff "%s:%d" env count
    ) fp.feat_counts;
  Buffer.contents buff

let debug_fp_str_log fp_str =
  let name, fp = S.split fp_str ~by:"\t" in
  Log.debug "name: %s" name;
  let toks = S.split_on_char ';' fp in
  (* sort atom envs. from least to most precise *)
  let sorted =
    L.sort (fun (s1: string) s2 ->
        let (r1: int) = S.count_char s1 ',' in
        let r2 = S.count_char s2 ',' in
        (* less specific atom envs come first *)
        if r1 < r2 then
          -1
        else if r1 > r2 then
          1
        else (* if same radius: alphabetical order *)
          compare s1 s2
      ) toks in
  L.iter (Log.debug "%s") sorted

(* unfolded counted RFP encoding *)
let encode_smiles_line max_radius line =
  let smi, name = BatString.split ~by:"\t" line in
  let mol_noH = Rdkit.__init__ ~smi () in
  (* this fingerprint needs all hydrogens *)
  let mol = Rdkit.add_hydrogens mol_noH () in
  let num_atoms = Rdkit.get_num_atoms mol () in
  let elements = Rdkit.get_elements mol () in
  let indexes = A.init num_atoms (fun i -> i) in
  { name;
    feat_counts =
      A.fold (fun acc0 a_i ->
          (* current atom's environments *)
          let buff = Buffer.create 128 in
          (* encode each atom using all radii; from 0 to max radius
             for this atom (furthest neighbor on the molecular graph) *)
          let dists = Rdkit.get_distances mol ~i:a_i () in          
          let radii =
            let r_max = 1 + (min max_radius (A.max dists)) in
            A.init r_max (fun i -> i) in
          A.fold (fun acc1 radius ->
              let atom_env = get_atom_env dists radius indexes elements in
              (* separate layers *)
              (if Buffer.length buff > 0 then
                 Buffer.add_char buff ','
              );
              (* chemical formula for this layer *)
              SMap.iter (fun k v ->
                  if v > 1 then
                    Printf.bprintf buff "%s%d" k v
                  else (* as in chemical formulae *)
                    Printf.bprintf buff "%s" k
                ) atom_env;
              let curr_env = Buffer.contents buff in
              let count = SMap.find_default 0 curr_env acc1 in
              SMap.add curr_env (count + 1) acc1
            ) acc0 radii
        ) SMap.empty indexes
  }

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
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let line = input_line input in
          incr reads;
          let fp = encode_smiles_line max_radius line in
          let fp_str = fp_to_string fp in
          (if verbose then
             debug_fp_str_log fp_str
           else
             Printf.fprintf output "%s\n" fp_str
          );
          incr writes
        done
      with End_of_file -> ()
    );
  let dt = Unix.gettimeofday() -. start in
  Log.info "wrote %#d molecules (%.2f Hz)" !writes ((float !writes) /. dt)

let () = main ()
