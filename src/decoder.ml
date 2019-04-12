
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
              --r-mode: create data and labels files for R\n\
              --max-feat <int>: maximum possible feature index\
                                (needed if --r-mode)\n"
       Sys.argv.(0);
     exit 1);
  let db_fn = CLI.get_string ["-i"] args in
  let r_output_mode = CLI.get_set_bool ["--r-mode"] args in
  let max_feat =
    if r_output_mode then CLI.get_int ["--max-feat"] args
    else -1 in
  let r_data_fn, r_labels_fn =
    if r_output_mode then
      let out_prfx = Filename.remove_extension db_fn in
      let data_fn, labels_fn =
        (out_prfx ^ "_data.csv", out_prfx ^ "_labels.csv") in
      Log.info "creating %s" data_fn;
      Log.info "creating %s" labels_fn;
      (data_fn, labels_fn)
    else ("/dev/null", "/dev/null") in
  CLI.finalize ();
  let all_lines = Utls.lines_of_file db_fn in
  match all_lines with
  | [] -> assert(false)
  | db_rad :: rest ->
    let all_mols = MSE_mol.of_lines rest in
    printf "%s\n" db_rad;
    let features_seen = Ht.create 11 in
    Utls.with_out_file r_data_fn (fun data_out ->
        Utls.with_out_file r_labels_fn (fun labels_out ->
            L.iteri (fun i mol ->
                let name = MSE_mol.get_name mol in
                let map = MSE_mol.get_map mol in
                (* feature values _MUST_ be printed out in increasing order
                   of feature ids; hence the IntMap we create *)
                let feat_counts =
                  StringMap.fold (fun feat count acc ->
                      let feat_id =
                        Ht.(find_default
                              features_seen feat (length features_seen)) in
                      Ht.replace features_seen feat feat_id;
                      IntMap.add feat_id count acc
                    ) map IntMap.empty in
                if r_output_mode then
                  begin
                    let line = csv_line_of_int_map max_feat feat_counts in
                    fprintf data_out "%s\n" line;
                    let label_int =
                      if String.starts_with name "active" then 1 else -1 in
                    fprintf labels_out (if i > 0 then "\t%d" else "%d") label_int
                  end
                else
                  printf "%s,0.0,[%s]\n" name (mop2d_line_of_int_map feat_counts)
              ) all_mols;
            fprintf labels_out "\n"
          )
      );
    Log.info "read %d molecules from %s (%d features)"
      (L.length all_mols) db_fn (Ht.length features_seen)

let () = main ()
