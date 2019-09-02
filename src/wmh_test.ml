
open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Utls = Molenc.Utls
module WMH = Molenc.WMH

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i encoded_molecules.txt\n" Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  (* read all molecules *)
  let molecules = FpMol.molecules_of_file input_fn in
  let nb_features = L.max (L.map FpMol.nb_features molecules) in
  let sparse_fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let bounds = WMH.bounds nb_features sparse_fingerprints in
  let idx2feat = WMH.lookup_table bounds in
  let rand_bound = A.length idx2feat in
  let feat2acc_bound = WMH.acc_bounds_table bounds in
  let dense_fingerprints = A.map (WMH.to_dense nb_features) sparse_fingerprints in
  let n = A.length sparse_fingerprints in
  Log.info "read %d molecules" n;
  (* compute Tani for many pairs (and compute scoring rate) *)
  Random.init 12345; (* seed PRNG *)
  let pairs = A.init 10_000 (fun _i -> (Random.int n, Random.int n)) in
  let dists = A.make 10_000 0.0 in
  let dt1, () = Utls.time_it (fun () ->
      A.iteri (fun i (i1, i2) ->
          let tani =
            Fp.tanimoto
              (A.unsafe_get sparse_fingerprints i1)
              (A.unsafe_get sparse_fingerprints i2) in
          A.unsafe_set dists i tani
        ) pairs
    ) in
  let tani_rate = (float n) /. dt1 in
  Log.info "Tani-rate: %.2f" tani_rate;
  let ks = [10; 20; 30; 40; 50; 100; 200; 500] in
  (* test the correctness and bench hashing and scoring speeds
     as a function of k (the number of hashes) *)
  L.iter (fun k ->
      let data_fn = sprintf "k_%03d.data" k in
      Utls.with_out_file data_fn (fun out ->
          (* hash them (and compute hashing rate) *)
          let seeds = WMH.get_seeds k in
          let rands = WMH.gen_rands seeds rand_bound in
          let dt0, hashes = Utls.time_it (fun () ->
              A.map (WMH.hash rands idx2feat feat2acc_bound) dense_fingerprints
            ) in
          Log.info "k: %d hashing-rate: %11.2f" k (float n /. dt0);
          (* compute estimated tani for the same pairs
             (and compute scoring rate) *)
          let est_dists = A.make 10_000 0.0 in
          let dt2, () = Utls.time_it (fun () ->
              A.iteri (fun i (i1, i2) ->
                  let m1 = A.unsafe_get hashes i1 in
                  let m2 = A.unsafe_get hashes i2 in
                  let tani = WMH.estimate_jaccard m1 m2 in
                  A.unsafe_set est_dists i tani
                ) pairs
            ) in
          let est_tani_rate = (float n) /. dt2 in
          (if est_tani_rate <= tani_rate
           then Log.warn
           else Log.info) "k: %d est-Tani-rate: %.2f accel: %.2f"
            k est_tani_rate (est_tani_rate /. tani_rate);
          A.iteri (fun i exact_dist ->
              let abs_error = abs_float (exact_dist -. est_dists.(i)) in
              fprintf out "%f %f %f\n" exact_dist est_dists.(i) abs_error
            ) dists;
          (* output maximum Tani error *)
          let diffs =
            A.map2 (fun d1 d2 -> abs_float (d1 -. d2)) dists est_dists in
          let max_error = A.max diffs in
          let avg_error = A.favg diffs in
          let med_error = Utls.list_medianf (A.to_list diffs) in
          Log.info "k: %d error(max, avg, med): %.2f %.2f %.2f"
            k max_error avg_error med_error
        )
    ) ks

let () = main ()
