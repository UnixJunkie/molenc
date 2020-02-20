
open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Log = Dolog.Log
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
  let ks = [40] in
  (* bench hashing and scoring speeds *)
  L.iter (fun k ->
      (* hash them (and compute hashing rate) *)
      let seeds = WMH.get_seeds k in
      let rands = WMH.gen_rands seeds rand_bound in
      let dt0, hashes = Utls.time_it (fun () ->
          A.map (WMH.hash rands idx2feat feat2acc_bound) dense_fingerprints
        ) in
      Log.info "k: %d hashing-rate: %.2f" k (float n /. dt0);
      (* compute estimated tani for the same pairs (and compute scoring rate) *)
      let dt2, _est_dists = Utls.time_it (fun () ->
          let res = A.make n 0.0 in
          for i = 0 to n - 1 do
            let i1 = Random.int n in
            let i2 = Random.int n in
            let m1 = A.get hashes i1 in
            let m2 = A.get hashes i2 in
            let tani = WMH.estimate_jaccard m1 m2 in
            A.set res i tani
          done;
          res) in
      let est_tani_rate = (float n) /. dt2 in
      Log.info "k: %d est-Tani-rate: %.2f" k est_tani_rate
    ) ks

let () = main ()
