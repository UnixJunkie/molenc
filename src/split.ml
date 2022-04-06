(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Split large files with many molecules into chunks of given size.
   .smi, .sdf and .mol2 are supported. *)

open Printf

module A = Array
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

type file_format = SMILES
                 | SDF
                 | Mol2

let file_format_of_filename fn =
  if S.ends_with fn ".smi" || S.ends_with fn ".smiles" then
    SMILES
  else if S.ends_with fn ".sdf" then
    SDF
  else if S.ends_with fn ".mol2" then
    Mol2
  else
    failwith
      ("Split.file_format_of_filename: not {.smi|.smiles|.sdf|.mol2}: %s" ^ fn)

exception Read_one

let read_one ff input =
  match ff with
  | Mol2 -> failwith "not implemented yet"
  | SMILES -> [input_line input]
  | SDF ->
    let res = ref [] in
    try
      while true do
        let line = input_line input in
        if line = "$$$$" then
          begin
            res := line :: !res;
            raise Read_one
          end
        else
          res := line :: !res;
      done;
      assert(false)
    with Read_one -> L.rev !res

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
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let _csize = CLI.get_int_def ["-c"] args 50 in
  CLI.finalize ();
  let _ff = file_format_of_filename input_fn in
  failwith "not implemented yet"

let () = main ()
