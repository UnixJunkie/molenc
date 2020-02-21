
(* initial tests of random coordinates sub-sampling wihout replacement
   (for fast but approximate Jaccard computation) *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -p <drop_frac:FLOAT> -i <FILE> [-n <repeats:INT>]\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let drop_p = CLI.get_float ["-p"] args in
  assert(drop_p > 0.0 && drop_p < 1.0);
  let nb_iter = CLI.get_int_def ["-n"] args 10_000 in
  CLI.finalize ();
  let alpha = 1.0 /. (1.0 -. drop_p) in
  (* read all molecules *)
  let molecules = FpMol.molecules_of_file input_fn in
  let nb_mols = L.length molecules in
  Log.info "nb_mols: %d" nb_mols;
  let fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let nb_features = L.max (L.map FpMol.nb_features molecules) in
  Log.info "nb_features: %d" nb_features;
  let feat_id_max = nb_features - 1 in
  let rand_feat_ids =
    let all_features = L.range 0 `To feat_id_max in
    L.shuffle all_features in
  let truncated =
    let n = Utls.ceili (drop_p *. (float nb_features)) in
    let to_drop = Ht.create n in
    let candidates = L.take n rand_feat_ids in
    L.iter (fun i ->
        Ht.add to_drop i ()
      ) candidates;
    A.map (Fp.drop_features to_drop) fingerprints in
  for _ = 1 to nb_iter do
    let i = Random.int nb_mols in
    let j = Random.int nb_mols in
    let fp_i = fingerprints.(i) in
    let fp_j = fingerprints.(j) in
    let sum_min, sum_max = Fp.sum_min_max fp_i fp_j in
    let exact_tani = Fp.tanimoto fp_i fp_j in
    let tfp_i = truncated.(i) in
    let tfp_j = truncated.(j) in
    let est_tani = Fp.tanimoto tfp_i tfp_j in
    let est_sum_min, est_sum_max = Fp.sum_min_max tfp_i tfp_j in
    let x = float sum_min /. float est_sum_min in
    let y = float sum_max /. float est_sum_max in
    printf "Tani:\t%.3f\t%.3f\t%.3f" exact_tani est_tani
      (abs_float (exact_tani -. est_tani));
    printf "\tn: %d\t%d\t%.3f\t%.3f\t%.3f" sum_min est_sum_min
      x alpha ((x -. alpha) /. x);
    printf "\tu: %d\t%d\t%.3f\t%.3f\t%.3f\n" sum_max est_sum_max
      y alpha ((y -. alpha) /. y)
  done

(* FBR: on stderr, print the average absolute error in each case *)

(* FBR: write unit test *)

let () = main ()
