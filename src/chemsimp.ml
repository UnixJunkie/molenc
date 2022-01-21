(* try to determine if the chemical similarity principle holds
   for a given (already-encoded) molecular dataset *)

open Printf

module A = BatArray
(* module BST = Bst.Bisec_tree.Make (Molenc.FpMol) *)
module CLI = Minicli.CLI
module Ht = BatHashtbl
module Log = Dolog.Log
module Utls = Molenc.Utls
module FpMol = Molenc.FpMol
module L = BatList

(* molecular distances are already normalized (Tanimoto dists);
   hence for visualization it is convenient if the IC50s are also normalized *)
let normalize mini maxi pIC50 =
  (pIC50 -. mini) /. (maxi -. mini)

let nearest_neighbor ht i n_mols =
  let nearest_d = ref 1.0 in
  let nearest = ref (-1) in
  for j = 0 to i - 1 do
    let d = Ht.find ht (j, i) in
    if d < !nearest_d then
      begin
        nearest_d := d;
        nearest := j
      end
  done;
  for j = i + 1 to n_mols - 1 do
    let d = Ht.find ht (i, j) in
    if d < !nearest_d then
      begin
        nearest_d := d;
        nearest := j
      end
  done;
  (!nearest, !nearest_d)

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
  let n_mols = A.length molecules_a in
  Utls.enforce (n_mols >= 50)
    (sprintf "<50 molecules in %s: %d" input_fn n_mols);
  (* constructing the bst would take too much time in that case *)
  Utls.enforce (n_mols <= 50_000)
    (sprintf ">50_000 molecules in %s: %d" input_fn n_mols);
  (* Log.info "indexing molecules...";
   * let bst = BST.(create 1 Two_bands molecules_a) in *)
  let min_pIC50, max_pIC50 = A.min_max (A.map FpMol.get_value molecules_a) in
  Log.info "pIC50_min/max: %g %g" min_pIC50 max_pIC50;
  let pairs = ref [] in
  (* let rng = BatRandom.State.make_self_init () in
   * for _i = 1 to num_draws do
   *   let m1 = molecules_a.(BatRandom.State.int rng n_mols) in
   *   let m2 = molecules_a.(BatRandom.State.int rng n_mols) in
   *   let dist = FpMol.dist m1 m2 in
   *   let act_abs_diff =
   *     let act1 = normalize min_pIC50 max_pIC50 (FpMol.get_value m1) in
   *     let act2 = normalize min_pIC50 max_pIC50 (FpMol.get_value m2) in
   *     abs_float (act1 -. act2) in
   *   Log.debug "%f %f\n" dist act_abs_diff;
   *   pairs := (dist, act_abs_diff) :: !pairs
   * done; *)
  (* initialize half of the Gram matrix *)
  let ht = Ht.create (n_mols * n_mols) in
  for i = 0 to n_mols - 2 do
    let m1 = molecules_a.(i) in
    for j = i + 1 to n_mols - 1 do
      let m2 = molecules_a.(j) in
      let d = FpMol.dist m1 m2 in
      Ht.add ht (i, j) d;
    done
  done;
  (* find the nearest of each; output distance and abs_act_diff *)
  for i = 0 to n_mols - 1 do
    let m1 = molecules_a.(i) in
    let m2_i, d = nearest_neighbor ht i n_mols in
    let m2 = molecules_a.(m2_i) in
    let act_abs_diff =
      let act1 = normalize min_pIC50 max_pIC50 (FpMol.get_value m1) in
      let act2 = normalize min_pIC50 max_pIC50 (FpMol.get_value m2) in
      abs_float (act1 -. act2) in
    Log.debug "%f %f\n" d act_abs_diff;
    pairs := (d, act_abs_diff) :: !pairs
  done;
  let dists, act_abs_diffs = L.split !pairs in
  let rho, t = Utls.spearman_l dists act_abs_diffs in
  Log.info "fn: %s N: %d rho: %.2f (p-value: %f)\n" input_fn num_draws rho t

let () = main ()
