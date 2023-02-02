(* Copyright (C) 2022, Francois Berenger

   Tsuda Laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Bounding Box Applicability Domain
   - support .AP files
   - support .csv files (output of molenc_lizard.py) *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IMap = BatMap.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

let is_AP_file fn =
  S.ends_with fn ".AP"

let is_csv_file fn =
  S.ends_with fn ".csv"

type mode = AP_files
          | CSV_files

let parse_AP_line line =
  Scanf.sscanf line "%s@[%s@]" (fun _name_ic50 feat_count_str ->
      let feat_counts = S.split_on_char ';' feat_count_str in
      L.rev_map (fun feat_count ->
          let feat, count = S.split ~by:":" feat_count in
          (int_of_string feat, int_of_string count)
        ) feat_counts
    )

(* array of molecular descriptors values *)
let parse_CSV_line line =
  (* the first feature is molecular name; let's skip this one *)
  let feature_strings = L.tl (S.split_on_char ',' line) in
  let strings = A.of_list feature_strings in
  A.map float_of_string strings

let ad_from_AP_file fn =
  LO.fold fn (fun acc1 line ->
      let feat_counts = parse_AP_line line in
      L.fold_left (fun acc2 (feat, count) ->
          let prev_max_count = IMap.find_default 0 feat acc2 in
          IMap.add feat (max prev_max_count count) acc2
        ) acc1 feat_counts
    ) IMap.empty

let ad_from_CSV_file fn =
  let nfields = ref 0 in
  let min_maxs = ref (Array.make 1 (0.0,0.0)) in
  LO.iteri fn (fun i line ->
      if S.starts_with line "#" then (* header line *)
        begin
          assert(i = 0); (* comment lines not allowed *)
          nfields := S.count_char line ','; (* ignore name (1st field) *)
          min_maxs := Array.make !nfields (infinity, neg_infinity);
        end
      else
        let features = parse_CSV_line line in
        assert(A.length features = !nfields);
        A.iteri (fun j x ->
            let curr_min, curr_max = !min_maxs.(j) in
            !min_maxs.(j) <- (min curr_min x, max curr_max x)
          ) features
    );
  !min_maxs

(* all feature counts are inside the bouding box *)
let is_in_AP_AD ad feat_counts =
  L.for_all (fun (feat, count) ->
      let max_count = IMap.find_default 0 feat ad in
      count <= max_count
    ) feat_counts

let apply_AP_AD_to_file ad fn =
  LO.filter fn (fun line ->
      let feat_counts = parse_AP_line line in
      is_in_AP_AD ad feat_counts
    )

let apply_CSV_AD_to_file ad fn =
  LO.filter fn (fun line ->
      if S.starts_with line "#" then
        true
      else
        let features = parse_CSV_line line in
        A.for_all2 (fun x (mini, maxi) ->
            x >= mini && x <= maxi
          ) features ad
    )

let string_from_AP_AD ad =
  let buff = Buffer.create 1024 in
  IMap.iter (fun k v ->
      Printf.bprintf buff "%d:%d " k v
    ) ad;
  Buffer.contents buff

let string_from_CSV_AD ad =
  let buff = Buffer.create 1024 in
  A.iter (fun (mini, maxi) ->
      Printf.bprintf buff "(%g,%g) " mini maxi
    ) ad;
  Buffer.contents buff

let main () =
  let _start = Unix.gettimeofday () in
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-tr|--train <train.{AP|csv}]: training set\n  \
              [-te|--test <test.{AP|csv}]: test set\n  \
              [-o output_fn]: output file (passed AD)\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let train_fn = CLI.get_string ["-tr";"--train"] args in
  let test_fn = CLI.get_string ["-te";"--test"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let file_type =
    match (is_AP_file train_fn, is_AP_file test_fn) with
    | true, true -> AP_files
    | _, _ ->
      match (is_csv_file train_fn, is_csv_file test_fn) with
      | true, true -> CSV_files
      | _, _ ->
        failwith "BBAD: only two .AP files or two .csv files allowed" in
  match file_type with
  | AP_files ->
    begin
      let ad = ad_from_AP_file train_fn in
      Log.debug "AD:";
      Log.debug "%s" (string_from_AP_AD ad);
      let before = LO.length test_fn in
      let in_AD = apply_AP_AD_to_file ad test_fn in
      let after = L.length in_AD in
      Log.info "before/after: %d/%d" before after;
      LO.lines_to_file output_fn in_AD
    end
  | CSV_files ->
    begin
      let ad = ad_from_CSV_file train_fn in
      Log.debug "AD:";
      Log.debug "%s" (string_from_CSV_AD ad);
      let before = LO.length test_fn in
      let in_AD = apply_CSV_AD_to_file ad test_fn in
      let after = L.length in_AD in
      Log.info "before/after: %d/%d" before after;
      LO.lines_to_file output_fn in_AD
    end

let () = main ()
