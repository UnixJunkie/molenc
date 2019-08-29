
open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Utls = Molenc.Utls
module WMH = Molenc.WMH
module Sfp = Molenc.Short_fingerprint

let main () =
  Log.color_on ();
  Log.set_log_level Log.DEBUG;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i encoded_molecules.txt -n nbits\n" Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let nbits = CLI.get_int ["-n"] args in
  (* read all molecules *)
  let molecules = FpMol.molecules_of_file input_fn in
  let nb_features = 1 + (L.max (L.map FpMol.max_feature_index molecules)) in
  let sparse_fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let bounds = WMH.bounds nb_features sparse_fingerprints in
  printf "bounds:\n";
  A.iteri (printf "%d:%d ") bounds;
  printf "\nEND\n";
  let idx2feat = WMH.lookup_table bounds in
  printf "idx2feat:\n";
  A.iter (printf "%d") idx2feat;
  printf "\nEND\n";
  let feat2acc_bound = WMH.acc_bounds_table bounds in
  printf "acc bounds:\n";
  A.iteri (printf "%d:%d ") feat2acc_bound;
  printf "\nEND\n";
  let dense_fingerprints = A.map (WMH.to_dense nb_features) sparse_fingerprints in
  A.iter (fun fp ->
      printf "%s\n" (WMH.string_of_dense fp)
    ) dense_fingerprints;
  let n = A.length sparse_fingerprints in
  Log.info "read %d molecules" n;
  let short_fps =
    A.map (Sfp.of_dense nbits idx2feat feat2acc_bound) dense_fingerprints in
  A.iter (fun fp ->
      printf "%s\n" (Sfp.to_string fp)
    ) short_fps

let () = main ()
