
(* Weighted Minwise Hashing in (amortized) Constant Time

   Shrivastava, A. (2016).
   Simple and efficient weighted minwise hashing.
   In Advances in Neural Information Processing Systems (pp. 1498-1506). *)

module A = BatArray
module BA = Bigarray
module BA1 = BA.Array1
module Fp = Fingerprint
module L = BatList

type dense = (int, BA.int8_unsigned_elt, BA.c_layout) BA1.t

(* FBR: later on, try using a BA for the hashed FP too *)
type hashed = int array

(* FBR: compute [k] in as a function of the precision we want
        (there is a formula in some papers) *)

let get_seeds k =
  let global_seed = 12345 in
  let rng = Random.State.make [|global_seed|] in
  let bound = (BatInt.pow 2 30) - 1 in
  Array.init k (fun _ -> Random.State.int rng bound)

(* FBR: put back unsafe where possible in whole file *)

(* convert the sparse Fp.t type into a dense array of small positive ints *)
let to_dense (feat_id_bound: int) (fp: Fp.t): dense =
  let res = BA1.create BA.int8_unsigned BA.C_layout feat_id_bound in
  BA1.fill res 0;
  let n = BA1.dim fp in
  let i = ref 0 in
  while !i < n do
    let k = BA1.get fp !i in
    let v = BA1.get fp (!i + 1) in
    assert(k < feat_id_bound);
    BA1.set res k v;
    i := !i + 2
  done;
  res

(* read the sparse fingerprints, update feat. val. bounds
 * if necessary *)
let update_bounds (bounds: int array) (fp: Fp.t): unit =
  let n = BA1.dim fp in
  let i = ref 0 in
  while !i < n do
    let k = BA1.get fp !i in
    let v = BA1.get fp (!i + 1) in
    bounds.(k) <- max (bounds.(k)) v;
    i := !i + 2
  done

(* compute the max value for each feature. *)
(* I.e. the columns' maximum if we put observations as rows
 * and features as columns in a data matrix *)
let bounds (max_feat_id: int) (train: Fp.t array): int array =
  let bounds = A.make max_feat_id 0 in
  A.iter (update_bounds bounds) train;
  bounds

(* create a lookup table from the bounds (max feature values) so that we
   can draw a single rand but still know which feature id. it corresponds to *)
let lookup_table (bounds: int array): int array =
  let total = A.sum bounds in
  let res = A.create total 0 in
  let j = ref 0 in
  A.iteri (fun i bound ->
      for _ = 1 to bound do
        res.(!j) <- i;
        incr j
      done
    ) bounds;
  res

let acc_bounds_table (bounds: int array): int array =
  let n = A.length bounds in
  let res = A.create n 0 in
  let acc = ref 0 in
  A.iteri (fun i bound ->
      res.(i) <- !acc;
      acc := !acc + bound
    ) bounds;
  res

(* in the paper, he defines is_green; but he samples until is_green becomes
 * true. It is more natural to sample while is_red *)
let is_red (arr: dense) (test_feat_id: int) (test_feat_val: int): bool =
  let feat_val = BA1.get arr test_feat_id in
  (feat_val = 0) || (test_feat_val > feat_val)

(* compute k hashes *)
let hash seeds idx2feat feat2acc_bound (dense_fp: dense): hashed =
  let k = A.length seeds in
  let rand_bound = A.length idx2feat in
  let res = A.make k 0 in
  for i = 0 to k - 1 do
    let misses = ref 0 in
    let seed = A.get seeds i in
    let rng = Random.State.make [|seed|] in
    let rand' = Random.State.int rng rand_bound in
    let test_feat_id = ref (idx2feat.(rand')) in
    let test_feat_val = ref (rand' - feat2acc_bound.(!test_feat_id)) in
    while is_red dense_fp !test_feat_id !test_feat_val do
      incr misses; (* in the paper: Hashes[i]++ *)
      let rand = Random.State.int rng rand_bound in
      test_feat_id := idx2feat.(rand);
      test_feat_val := rand - feat2acc_bound.(!test_feat_id)
    done;
    res.(i) <- !misses
  done;
  res

let estimate_jaccard (hash1: hashed) (hash2: hashed): float =
  let res = ref 0 in
  let k = A.length hash1 in
  assert(A.length hash2 = k);
  for i = 0 to k - 1 do
    if hash1.(i) = hash2.(i) then
      incr res
  done;
  (float !res) /. (float k)

let estimate_distance (h1: hashed) (h2: hashed): float =
  1.0 -. (estimate_jaccard h1 h2)
