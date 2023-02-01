(* Copyright (C) 2022, Francois Berenger

   Tsuda Laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Bounding Box Applicability Domain
   - support .AP files
   - support .csv files (output of molenc_lizard.py) *)

open Printf

module CLI = Minicli.CLI
module Ht = BatHashtbl
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

let main () =
  let _start = Unix.gettimeofday () in
  Log.(set_prefix_builder short_prefix_builder);
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename.ph4>]: input file\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let train_fn = CLI.get_string ["-tr";"--train"] args in
  let test_fn = CLI.get_string ["-te";"--test"] args in
  let file_type =
    match (is_AP_file train_fn, is_AP_file test_fn) with
    | true, true -> AP_files
    | _, _ ->
      match (is_csv_file train_fn, is_csv_file test_fn) with
      | true, true -> CSV_files
      | _, _ ->
        failwith "BBAD: only two .AP files or two .csv files allowed" in
  match file_type with
  | AP_files -> failwith "not implemented yet"
  | CSV_files -> failwith "not implemented yet"

let () = main ()
