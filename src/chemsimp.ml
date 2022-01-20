(* try to determine if the chemical similarity principle holds
   for a given (already-encoded) molecular dataset *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Log = Dolog.Log
module Utls = Molenc.Utls
module FpMol = Molenc.FpMol
module L = BatList

(* molecular distances are already normalized (Tanimoto dists);
   hence for visualization it is convenient if the IC50s are also normalized *)
let normalize mini maxi pIC50 =
  (pIC50 -. mini) /. (maxi -. mini)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file (.AP format)\n  \
               [-n <int>]: number of pairs to draw; default=1000\n  \
               [-v]: verbose/debug mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let num_draws = CLI.get_int_def ["-n"] args 1000 in
  let verbose = CLI.get_set_bool ["-v"] args in
  CLI.finalize(); (* ------------------------------------------------------- *)
  if verbose then
    Log.(set_log_level DEBUG)
  ;
  let molecules_a = A.of_list (FpMol.molecules_of_file input_fn) in
  let min_pIC50, max_pIC50 = A.min_max (A.map FpMol.get_value molecules_a) in
  Log.info "pIC50_min/max: %g %g" min_pIC50 max_pIC50;
  let n_mols = A.length molecules_a in
  Utls.enforce (n_mols >= 50)
    (sprintf "less than 50 molecules in %s: %d" input_fn n_mols);
  let rng = BatRandom.State.make_self_init () in
  let pairs = ref [] in
  for _i = 1 to num_draws do
    let m1 = molecules_a.(BatRandom.State.int rng n_mols) in
    let m2 = molecules_a.(BatRandom.State.int rng n_mols) in
    let dist = FpMol.dist m1 m2 in
    let act_abs_diff =
      let act1 = normalize min_pIC50 max_pIC50 (FpMol.get_value m1) in
      let act2 = normalize min_pIC50 max_pIC50 (FpMol.get_value m2) in
      abs_float (act1 -. act2) in
    Log.debug "%f %f\n" dist act_abs_diff;
    pairs := (dist, act_abs_diff) :: !pairs
  done;
  let dists, act_abs_diffs = L.split !pairs in
  let rho, t = Utls.spearman_l dists act_abs_diffs in
  Log.info "N: %d rho: %.2f (p-value: %f)\n" num_draws rho t

let () = main ()
