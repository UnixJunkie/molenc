
(* molecular decoder
   convert the MSE txt format to a 1mop2d-compatible text format
   optionally, this can also decode towards an R friendly CSV format *)

open Printf

module CLI = Minicli.CLI
module L = MyList
module Ht = BatHashtbl
module String = BatString
module StringMap = BatMap.String
module IntMap = BatMap.Int

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

let iwn_line_of_int_map map =
  let buff = Buffer.create 11 in
  let start = ref true in
  let total = float (IntMap.fold (fun _k v acc -> acc + v) map 0) in
  IntMap.iter (fun k v ->
      let scaled = (float v) /. total in
      if !start then
        (bprintf buff "%d:%f" k scaled;
         start := false)
      else
        bprintf buff " %d:%f" k scaled
    ) map;
  Buffer.contents buff

let csv_line_of_int_map max_feat map =
  let buff = Buffer.create 11 in
  let max_key, _v = IntMap.max_binding map in
  if max_key > max_feat then
    failwith (sprintf
                "Decoder.csv_line_of_int_map: max_key (%d) > max_feat (%d)"
                max_key max_feat);
  for k = 0 to max_feat do
    let v = IntMap.find_default 0 k map in
    bprintf buff (if k > 0 then "\t%d" else "%d") v
  done;
  Buffer.contents buff

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i db\n\
              -i <filename>: encoded molecules database\n\
              -o <filename>: where to write decoded molecules\n\
              --iwn: perform Instance-Wise Normalisation\n\
              --bitstring <filename>: output FPs as bitstrings\n"
       Sys.argv.(0);
     exit 1);
  let db_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_bin_fn = CLI.get_string_opt ["--bitstring"] args in
  let normalize = CLI.get_set_bool ["--iwn"] args in
  CLI.finalize ();
  let all_lines = Utls.lines_of_file db_fn in
  match all_lines with
  | [] -> assert(false)
  | db_rad :: rest ->
    let all_mols = MSE_mol.of_lines rest in
    let nb_mols = L.length all_mols in
    Log.info "%s" db_rad;
    let feat_to_id = Ht.create nb_mols in
    let feat_id_to_max_count = Ht.create nb_mols in
    let mol_name_idx_to_feat_counts = Ht.create nb_mols in
    Utls.with_out_file output_fn (fun out ->
        L.iteri (fun i mol ->
            let name = MSE_mol.get_name mol in
            let map = MSE_mol.get_map mol in
            (* feature values _MUST_ be printed out in increasing
                       order of feature ids; hence the IntMap we create *)
            let feat_counts =
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
                ) map IntMap.empty in
            Ht.add mol_name_idx_to_feat_counts (name, i) feat_counts;
            if normalize then
              let label =
                if BatString.starts_with name "active" then 1 else -1 in
              let line = iwn_line_of_int_map feat_counts in
              fprintf out "%+d %s\n" label line
            else
              let line = mop2d_line_of_int_map feat_counts in
              fprintf out "%s,0.0,[%s]\n" name line
          ) all_mols;
      );
    (* FBR: the dictionary should go in its own file; or at least as comment
     * lines in the output file *)
    (* output dictionary and max_counts *)
    let incr_feat_ids =
      let feat_ids' = Ht.to_list feat_to_id in
      L.sort (fun (_feat1, id1) (_feat2, id2) ->
          BatInt.compare id1 id2
        ) feat_ids' in
    printf "#DICTIONARY:\n";
    let max_bitwidth = ref 0 in
    let total_bits_required =
      L.fold_left (fun acc (feature, id) ->
          let max_count = Ht.find feat_id_to_max_count id in
          max_bitwidth := max !max_bitwidth max_count;
          printf "#%d:%d:%s\n" id max_count feature;
          acc + max_count
        ) 0 incr_feat_ids in
    let nb_features = Ht.length feat_to_id in
    Log.info "read %d molecules from %s" nb_mols db_fn;
    Log.info "#features: %d largest_field: %d #bits/ligand: %d"
      nb_features !max_bitwidth total_bits_required;
    (* output in binary strings format *)
    match maybe_bin_fn with
    | None -> ()
    | Some bin_fn ->
      Utls.with_out_file bin_fn (fun bin_out ->
          L.iteri (fun i mol ->
              (* molecule name *)
              let name = MSE_mol.get_name mol in
              fprintf bin_out "%s,0.0," name;
              let feat_counts =
                Ht.find mol_name_idx_to_feat_counts (name, i) in
              (* bitstring representation *)
              for feat_id = 0 to nb_features - 1 do
                let nb_bits = Ht.find feat_id_to_max_count feat_id in
                let count = IntMap.find_default 0 feat_id feat_counts in
                let ones = String.make count '1' in
                let zeroes = String.make (nb_bits - count) '0' in
                fprintf bin_out "%s%s" ones zeroes
              done;
              (* terminate line *)
              fprintf bin_out "\n"
            ) all_mols;
        )

let () = main ()
