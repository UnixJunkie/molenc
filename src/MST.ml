(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   Minimum Spanning Tree using the dataset's Gram matrix. *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename>: encode molecules file\n  \
              -o <filename>: output file\n"
       Sys.argv.(0);
     exit 1);
  let _input_fn = CLI.get_string ["-i"] args in
  CLI.finalize ();
  failwith "not implemented yet"

let () = main ()
