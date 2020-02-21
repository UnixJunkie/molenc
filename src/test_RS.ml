
(* initial tests of random coordinates sub-sampling wihout replacement
   (for fast but approximate Jaccard replacement) *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

(* FBR: make this a one shot program *)

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
  let nb_iter = CLI.get_int_def ["-n"] args 10_000 in
  CLI.finalize ();
  let alpha = 1.0 /. (1.0 -. drop_p) in
  (* read all molecules *)
  let molecules = FpMol.molecules_of_file input_fn in
  let nb_mols = L.length molecules in
  let fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let nb_features = L.max (L.map FpMol.nb_features molecules) in
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
  (* let abs_err_acc = ref 0.0 in *)
  printf "#drop_p: %f\n" drop_p;
  printf "#exact   est1     err1     est2     err2     alpha    Ds\n";
  (* "    0.205479 0.163265 0.042214 0.114943 0.090537 1.428571" *)
  for _ = 1 to nb_iter do
    let i = Random.int nb_mols in
    let j = Random.int nb_mols in
    let fp_i = fingerprints.(i) in
    let fp_j = fingerprints.(j) in
    let norm_fp_i = Fp.sum_values fp_i in
    let norm_fp_j = Fp.sum_values fp_j in
    let exact_tani = Fp.tanimoto fp_i fp_j in
    let tfp_i = truncated.(i) in
    let tfp_j = truncated.(j) in
    (* D_s in the paper *)
    let ds =
      (float nb_features) /.
      (float (min (Fp.nb_features tfp_i) (Fp.nb_features tfp_j))) in
    let fully_estimated_tani = Fp.tanimoto tfp_i tfp_j in
    let err_1 = abs_float (exact_tani -. fully_estimated_tani) in
    let est_sum_min, _est_sum_max = Fp.sum_min_max tfp_i tfp_j in
    let alpha_est_sum_min = alpha *. (float est_sum_min) in
    (* FBR: we should max bound this one to 1.0 *)
    let part_est_tani =
      alpha_est_sum_min /.
      ((float (norm_fp_i + norm_fp_j)) +. alpha_est_sum_min) in
    let err_2 = abs_float (exact_tani -. part_est_tani) in
    printf "%f %f %f %f %f %f %f\n"
      exact_tani
      fully_estimated_tani err_1
      part_est_tani err_2 alpha ds
  done

(* FBR: on stderr, print the average absolute error in each case *)

(* FBR: write unit test *)

let () = main ()
