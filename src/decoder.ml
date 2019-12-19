(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* molecular decoder
   convert the MSE txt format to a 1mop2d-compatible text format
   optionally, this can also decode towards an R friendly CSV format *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = BatHashtbl
module IntMap = BatMap.Int
module L = BatList
module Log = Dolog.Log
module MSE_mol = Molenc.MSE_mol
module Norm = Molenc.Norm
module String = BatString
module StringMap = BatMap.String
module Utls = Molenc.Utls

let mop2d_line_of_int_map map =
  let buff = Buffer.create 11 in
  let start = ref true in
  IntMap.iter (fun k v ->
      if !start then
        (bprintf buff "%d:%d" k v;
         start := false)
      else
        bprintf buff ";%d:%d" k v
    ) map;
  Buffer.contents buff

let iwn_line_of_int_map style map =
  let buff = Buffer.create 11 in
  let start = ref true in
  let total = Norm.map_norm style map in
  IntMap.iter (fun k' v ->
      let k = k' + 1 in
      let scaled = (float v) /. total in
      if !start then
        (bprintf buff "%d:%f" k scaled;
         start := false)
      else
        bprintf buff " %d:%f" k scaled
    ) map;
  Buffer.contents buff

(* reconstruct map/dico feat->featId from given file *)
let dictionary_from_file fn =
  let feat_to_id = Ht.create 10_000 in
  Utls.iter_on_lines_of_file fn (fun line ->
      if not (BatString.starts_with line "#") then
        Scanf.sscanf line "%d %d %s" (fun id _max_count feat ->
            (* the binding defined in the dictionary should be unique *)
            assert(not (Ht.mem feat_to_id feat));
            Ht.add feat_to_id feat id
          )
    );
  feat_to_id

type filename = string
type dico_mode = Read_from of filename
               | Write_to of filename

let feat_counts_dico_RW feat_to_id feat_id_to_max_count map =
  StringMap.fold (fun feat count acc ->
      let curr_nb_feats = Ht.length feat_to_id in
      let feat_id =
        Ht.find_default
          feat_to_id feat curr_nb_feats in
      Ht.replace feat_to_id feat feat_id;
      let prev_max_count =
        Ht.find_default feat_id_to_max_count feat_id 0 in
      Ht.replace feat_id_to_max_count
        feat_id (max prev_max_count count);
      IntMap.add feat_id count acc
    ) map IntMap.empty

let feat_counts_dico_RO feat_to_id map =
  StringMap.fold (fun feat count acc ->
      try
        let feat_id = Ht.find feat_to_id feat in
        IntMap.add feat_id count acc
      with Not_found -> acc (* ignored: feature not seen in training set
                               as a warning: we could count the number
                               of unknown features per molecule but... *)
    ) map IntMap.empty

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i db\n  \
              -i <filename>: encoded molecules database\n  \
              -o <filename>: where to write decoded molecules\n  \
              [--norm {l1|max}]: perform Instance-Wise Normalisation\n  \
              [-d <filename>]: read feature dictionary from file\n  \
              [-n <int>]: max number of parallel jobs\n"
       Sys.argv.(0);
     exit 1);
  let db_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-n"] args 1 in
  let dico = match CLI.get_string_opt ["-d"] args with
    | None -> Write_to (output_fn ^ ".dix")
    | Some fn -> Read_from fn in
  let maybe_norm =
    Utls.may_apply Norm.of_string (CLI.get_string_opt ["--norm"] args) in
  (* the normalized output format is for liblinear.
     Since molecule names are lost: it is FORBIDDEN to change
     the order between input and output files.
     Hence, parallelization is forbidden. *)
  Utls.enforce ((Option.is_none maybe_norm) || (nprocs = 1))
    "Decoder: normalized output only allowed in sequential mode";
  CLI.finalize ();
  Log.info "reading molecules...";
  let db_rad = Utls.get_first_line db_fn in
  let nb_mols =
    let cmd = sprintf "egrep -c '^#' %s" db_fn in
    (* rm header line from the count *)
    int_of_string (Utls.get_command_output cmd) - 1 in
  let feat_to_id = match dico with
    | Read_from dico_fn -> dictionary_from_file dico_fn
    | Write_to _fn -> begin
        Utls.enforce (nprocs = 1)
          "Decoder: writing the feature dictionary only works \
           in sequential mode";
        Ht.create nb_mols
      end in
  let feat_id_to_max_count = Ht.create nb_mols in
  let mol_count = ref 0 in
  Utls.with_infile_outfile db_fn output_fn (fun input output ->
      let header = input_line input in (* skip header comment line *)
      assert(BatString.starts_with header "#");
      try
        while true do
          let some_lines = MSE_mol.get_lines input in
          let mol = MSE_mol.read_one some_lines in
          let name = MSE_mol.get_name mol in
          let map = MSE_mol.get_map mol in
          (* feature values _MUST_ be printed out in increasing
             order of feature ids; hence the IntMap we create *)
          let feat_counts = match dico with
            | Write_to _ ->
              feat_counts_dico_RW feat_to_id feat_id_to_max_count map
            | Read_from _ ->
              feat_counts_dico_RO feat_to_id map in
          (match maybe_norm with
           | Some norm ->
             let label =
               if BatString.starts_with name "active" then 1 else -1 in
             let line = iwn_line_of_int_map norm feat_counts in
             fprintf output "%+d %s\n" label line
           | None ->
             let bitstring = mop2d_line_of_int_map feat_counts in
             fprintf output "%s,0.0,[%s]\n" name bitstring
          );
          if (!mol_count mod 1000) = 0 then
            printf "done: %d/%d\r%!" (!mol_count + 1) nb_mols;
          incr mol_count
        done
      with End_of_file ->
        Log.info "finished reading %s" db_fn
    );
  let incr_feat_ids =
    let feat_ids' = Ht.to_list feat_to_id in
    L.sort (fun (_feat1, id1) (_feat2, id2) ->
        BatInt.compare id1 id2
      ) feat_ids' in
  (* output dictionary and max_counts *)
  (match dico with
   | Read_from _ -> () (* read-only dico *)
   | Write_to dico_fn ->
     Utls.with_out_file dico_fn (fun out ->
         fprintf out "%s\n#featId maxCount featStr\n" db_rad;
         L.iter (fun (feature, id) ->
             let max_count = Ht.find feat_id_to_max_count id in
             fprintf out "%d %d %s\n" id max_count feature
           ) incr_feat_ids
       );
     Log.info "dictionary written to %s" dico_fn
  );
  let nb_features = Ht.length feat_to_id in
  Log.info "read %d molecules from %s" nb_mols db_fn;
  Log.info "#features: %d" nb_features

let () = main ()
