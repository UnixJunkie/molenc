(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * UHD encoder *)

(* FBR:TODO
   - compute dictionaries for ChEMBL-34: r=1, r=2, r=3 and r=+inf; DONE; process it
   - output number of times feature was seen in the dictionary file
     - this info is not read when loading in a dictionary; but it can be used
       to reorder features; from most frequent to least frequent
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

(* is the feature dictionary read in or written out *)
type dict_mode = Input of string
               | Output of string

type fp = { name: string;
            feat_counts: int SMap.t }

(* chemical formula at [radius] bonds away from center atom *)
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
  let elements = Rdkit.get_elements mol () in
  let atom_envs =
    (* gather all heavy atom environments *)
    A.fold_lefti (fun acc i elt ->
        if elt = "H" then
          acc (* only consider heavy atoms to reduce dimensionality of the FP *)
        else
          let dists = Rdkit.get_distances mol ~i () in
          (environments_for_atom dists fp_radius elements) :: acc
      ) [] elements in
  (* count them *)
  L.fold_left (fun acc0 l ->
      L.fold_left (fun acc1 formula ->
          let count = SMap.find_default 0 formula acc1 in
          SMap.add formula (count + 1) acc1
        ) acc0 l
    ) SMap.empty atom_envs

(* unfolded counted fingerprint encoding *)
let encode_smiles_line max_radius line =
  let smi, name = BatString.split ~by:"\t" line in
  let mol_noH = Rdkit.__init__ ~smi () in
  { name; feat_counts = environments_for_molecule max_radius mol_noH }

let fp_string_output mode out dict fp =
  fprintf out "%s,0.0,[" fp.name;
  let feat_counts = match mode with
    | Output _ -> (* writable dict *)
      (* feature index 0 is a reserved value for later users of the same dict:
         unkown feature *)
      let feat_i = ref (1 + Ht.length dict) in
      SMap.fold (fun feat count acc ->
          let feat_idx, seen =
            try Ht.find dict feat
            with Not_found ->
              let idx = !feat_i in
              Ht.add dict feat (idx, 0);
              incr feat_i;
              (idx, 0) in
          Ht.replace dict feat (feat_idx, seen + count);
          IMap.add feat_idx count acc
        ) fp.feat_counts IMap.empty
    | Input _ -> (* read-only dict *)
      SMap.fold (fun feat count acc ->
          try IMap.add (fst (Ht.find dict feat)) count acc
          with Not_found ->
            (* let () = Log.warn "unknown feat: %s" feat in *)
            let prev_count = IMap.find_default 0 0 acc in
            IMap.add 0 (1 + prev_count) acc
        ) fp.feat_counts IMap.empty in
  (* warn if unknown features *)
  (try Log.warn "%s: %d unknown features" fp.name (IMap.find 0 feat_counts)
   with Not_found -> ());
  let started = ref false in
  IMap.iter (fun feat count ->
      if !started then
        fprintf out ";%d:%d" feat count
      else
        (started := true;
         fprintf out "%d:%d" feat count)
    ) feat_counts;
  fprintf out "]\n"

let dico_format_header = "UHD-1.0.0"

let dico_from_file fn =
  Log.info "reading %s" fn;
  let lines = LO.lines_of_file fn in
  match lines with
  | [] -> (Log.fatal "Molenc_UHD.dico_from_file: %s is empty" fn;
           exit 1)
  | header :: body ->
    (* The Ultra High Dimensional fp? *)
    let () = assert(header = dico_format_header) in
    let dict = Ht.create (L.length body) in
    L.iter (fun line ->
        Scanf.sscanf line "%s@\t%d\t%d" (fun feat idx count ->
            Ht.add dict feat (idx, count)
          )
      ) body;
    Log.info "%d features" (Ht.length dict);
    dict

let dico_to_file dict fn =
  Log.info "writing feature dictionary to %s (%#d features)"
    fn (Ht.length dict);
  let kvs = Ht.to_list dict in
  (* sort by incr. feat index *)
  let kvs =
    L.sort (fun (_f1, (i, _count1)) (_f2, (j, _count2)) ->
        BatInt.compare i j) kvs in
  LO.with_out_file fn (fun output ->
      fprintf output "%s\n" dico_format_header;
      L.iter (fun (feat, (idx, count)) ->
          fprintf output "%s\t%d\t%d\n" feat idx count
        ) kvs
    )

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
               [-d <dict.dix>: use existing feature dictionary\n  \
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
  let mode = match CLI.get_string_opt ["-d"] args with
    | None -> Output (input_fn ^ ".dix")
    | Some dict_fn -> Input dict_fn in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 200 in
  let force = CLI.get_set_bool ["-f"] args in
  let max_radius = match CLI.get_int_opt ["-m"] args with
    | None -> max_int
    | Some x -> x in
  let verbose = CLI.get_set_bool ["-v"] args in
  CLI.finalize (); (* ------------------------------------------------------ *)
  (if verbose then
     Log.(set_log_level DEBUG)
  );
  let already_out_fn = Sys.file_exists output_fn in
  (if already_out_fn then
     (Log.warn "Molenc_UHD: output file exists: %s" output_fn;
      if not force then
        exit 1
     );
  );
  let writes = ref 0 in
  let dict = match mode with
    | Output _ -> Ht.create 20_011
    | Input fn -> dico_from_file fn in
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      Parany.run ~csize ~preserve:true nprocs
        ~demux:(fun () ->
            try input_line input
            with End_of_file -> raise Parany.End_of_input)
        ~work:(encode_smiles_line max_radius)
        ~mux:(function fp ->
            fp_string_output mode output dict fp;
            incr writes;
            if !writes mod 1_000 = 0 then
              eprintf "wrote: %#d\r%!" !writes)
    );
  let dt = Unix.gettimeofday() -. start in
  Log.info "wrote %#d molecules (%.2f Hz)" !writes ((float !writes) /. dt);
  match mode with
  | Output dict_fn -> dico_to_file dict dict_fn
  | Input _ -> ()

let () = main ()
