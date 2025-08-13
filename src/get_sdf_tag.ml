(* Efficiently extract tags from .sdf[.gz] files *)

open Printf

module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString
module StringSet = BatSet.String

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.{sdf[.gz]|mol2|smi|ph4} \
              {-names \"mol1,mol2,...\"|-f names_file} [-v]\n  \
              -i <filename>: molecules input file\n  \
              [-o <filename>]: molecules output file (default=stdout)\n  \
              [-names <string>,<string>,...]: molecule names\n  \
              [-f <filename>]: get molecule names from file\n  \
              [-if <filename>,<filename>,...]: several molecule input files\n  \
              [--force]: overwrite existing db file(s), if any\n  \
              [--no-index]: do not create db file(s)\n"
       Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  Log.set_log_level (if verbose then Log.DEBUG else Log.INFO);
  Log.color_on ();
  let _maybe_output_fn = CLI.get_string_opt ["-o"] args in
  ()
