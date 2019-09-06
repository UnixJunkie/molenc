(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Clustering using the Butina algorithm.
   -> output cluster centers
   -> output cluster membership for each molecule *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Utls = Molenc.Utls

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_list mols =
    create 1 Two_bands (Array.of_list mols)
end

let verbose = ref false

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: molecules to filter (\"database\")\n  \
               -o <filename>: output file\n  \
               [-t <float>]: Tanimoto threshold (default=1.0)\n  \
               [-e <filename>]: molecules to exclude from the ones to filter\n  \
               [-a <filename>]: molecules to annotate the ones to filter\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let _output_fn = CLI.get_string ["-o"] args in
  let threshold = CLI.get_float_def ["-t"] args 0.9 in
  assert(threshold >= 0.0 && threshold <= 1.0);
  verbose := CLI.get_set_bool ["-v"] args;
  CLI.finalize ();
  let mols_to_filter = FpMol.molecules_of_file input_fn in
  let total = L.length mols_to_filter in
  Log.info "read %d from %s" total input_fn

let () = main ()
