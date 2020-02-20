
open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls
module WMH = Molenc.WMH

let print_array title a =
  printf "%s:" title;
  for i = 0 to (A.length a) - 1 do
    printf " %d" a.(i)
  done;
  printf "\n"

let printi_array title a =
  printf "%s:" title;
  for i = 0 to (A.length a) - 1 do
    printf " %d:%d" i a.(i)
  done;
  printf "\n"

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  (* read all molecules *)
  let molecules =
    L.mapi FpMol.parse_one
      ["m0,0.0,[2:1;3:1]";
       "m1,0.0,[0:2;1:2;2:2;3:2]"] in
  let nb_features = L.max (L.map FpMol.nb_features molecules) in
  printf "nb_features: %d\n" nb_features;
  let sparse_fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let bounds = WMH.bounds nb_features sparse_fingerprints in
  print_array "bounds" bounds;
  let idx2feat = WMH.lookup_table bounds in
  printi_array "idx2feat" idx2feat;
  let rand_bound = A.length idx2feat in
  printf "bound: %d\n" rand_bound;
  let feat2acc_bound = WMH.acc_bounds_table bounds in
  print_array "feat2acc_bound" feat2acc_bound;
  let dense_fingerprints = A.map (WMH.to_dense nb_features) sparse_fingerprints in
  (* k = 1 *)
  let rands = [|0;1;2;3;5;7;4|] in (* only last one should hit *)
  let hash = WMH.hash [|rands|] idx2feat feat2acc_bound dense_fingerprints.(0) in
  let hash_val = hash.(0) in
  assert(hash_val = 6);
  let rands = [|0;1;2;3;7;5;6|] in (* only last one should hit *)
  let hash = WMH.hash [|rands|] idx2feat feat2acc_bound dense_fingerprints.(0) in
  let hash_val = hash.(0) in
  assert(hash_val = 6);
  let rands = [|4|] in (* 1st one should hit *)
  let hash = WMH.hash [|rands|] idx2feat feat2acc_bound dense_fingerprints.(0) in
  let hash_val = hash.(0) in
  assert(hash_val = 0);
  let rands = [|6|] in (* 1st one should hit *)
  let hash = WMH.hash [|rands|] idx2feat feat2acc_bound dense_fingerprints.(0) in
  let hash_val = hash.(0) in
  assert(hash_val = 0)

let () = main ()
