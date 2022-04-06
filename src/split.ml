(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo, Japan.

   Split large files with many molecules into chunks of given size. *)

open Printf

module A = Array
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename>]: SMILES input file\n  \
              -o <filename>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [--bbad <fn1[,fn2[,...]]>]: previously computed BBAD\n  \
              to apply as filter\n  \
              will compute the union of ADs if several filenames\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let _input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  CLI.finalize ();
  failwith "not implemented yet"

let () = main ()
