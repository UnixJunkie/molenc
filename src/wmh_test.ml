
open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Utls = Molenc.Utls
module WMH = Molenc.WMH

(* FBR: put a file of 100_000 encoded molecules in data/ *)

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
  let sparse_fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  (* 749: max feature index in test data file data/100k_mols_std_01.txt *)
  let dense_fingerprints = A.map (WMH.to_dense 749) sparse_fingerprints in
  let n = A.length sparse_fingerprints in
  Log.info "read %d molecules" n;
  (* compute Tani for many pairs (and compute scoring rate) *)
  Random.init 12345; (* seed PRNG *)
  let dt1, dists = Utls.time_it (fun () ->
      let res = A.make n 0.0 in
      for i = 0 to n - 1 do
        let i1 = Random.int n in
        let i2 = Random.int n in
        let m1 = A.get sparse_fingerprints i1 in
        let m2 = A.get sparse_fingerprints i2 in
        let tani = Fp.tanimoto m1 m2 in
        A.set res i tani
      done;
      res
    ) in
  Log.info "Tani-rate: %.3f" (float n /. dt1);
  let ks = [10; 20; 50; 100; 200; 500] in
  (* test the correctness and bench hashing and scoring speeds
     as a function of k (the number of hashes) *)
  L.iter (fun k ->
      (* hash them (and compute hashing rate) *)
      let seeds = WMH.get_seeds k in
      Gc.full_major ();
      let dt0, hashes = Utls.time_it (fun () ->
          A.map (WMH.hash seeds) dense_fingerprints
        ) in
      Log.info "k: %d hashing-rate: %.3f" k (float n /. dt0);
      (* compute estimated tani for the same pairs (and compute scoring rate) *)
      Gc.full_major ();
      Random.init 12345; (* seed PRNG *)
      let dt2, est_dists = Utls.time_it (fun () ->
          let res = A.make n 0.0 in
          for i = 0 to n - 1 do
            let i1 = Random.int n in
            let i2 = Random.int n in
            let m1 = A.get hashes i1 in
            let m2 = A.get hashes i2 in
            let tani = WMH.estimate_jaccard m1 m2 in
            A.set res i tani
          done;
          res
        ) in
      Log.info "k: %d est-Tani-rate: %.3f" k (float n /. dt2);
      (* output maximum Tani error *)
      let diffs = A.map2 (fun d1 d2 -> abs_float (d1 -. d2)) dists est_dists in
      let max_abs_error = A.max diffs in
      Log.info "k: %d max-abs-error: %f" k max_abs_error
    ) ks

let () = main ()
