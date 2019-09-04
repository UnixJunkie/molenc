
(* Weighted Minwise Hashing in (amortized) Constant Time

   Shrivastava, A. (2016).
   Simple and efficient weighted minwise hashing.
   In Advances in Neural Information Processing Systems (pp. 1498-1506). *)

open Printf

module A = BatArray
module BA = Bigarray
module BA1 = BA.Array1
module Fp = Fingerprint
module L = BatList

type dense = (int, BA.int8_unsigned_elt, BA.c_layout) BA1.t

type hashed = int array

let get_seeds k =
  let global_seed = 314159265 in
  let rng = Random.State.make [|global_seed|] in
  let bound = (BatInt.pow 2 30) - 1 in
  Array.init k (fun _ -> Random.State.int rng bound)

(* convert the sparse Fp.t type into a dense array of small positive ints *)
let to_dense (feat_id_bound: int) (fp: Fp.t): dense =
  let res = BA1.create BA.int8_unsigned BA.C_layout feat_id_bound in
  BA1.fill res 0;
  let n = BA1.dim fp in
  let i = ref 0 in
  while !i < n do
    let k = BA1.get fp !i in
    let v = BA1.get fp (!i + 1) in
    assert(k < feat_id_bound && v < 256);
    BA1.set res k v;
    i := !i + 2
  done;
  res

let string_of_dense x =
  let n = BA1.dim x in
  let buff = Buffer.create 80 in
  for i = 0 to n - 1 do
    let j = BA1.get x i in
    bprintf buff " %d:%d" i j
  done;
  Buffer.contents buff

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

(* (\* in the paper, he defines is_green; but he samples until is_green becomes
 *  * true. It is more natural to sample while is_red *\)
 * let is_red (arr: dense) (test_feat_id: int) (test_feat_val: int): bool =
 *   test_feat_val >= (BA1.unsafe_get arr test_feat_id) *)

(* pre-generate non-repeating random number sequences for later *)
let gen_rands seeds rand_bound =
  A.map (fun seed ->
      let rng = Random.State.make [|seed|] in
      (* all ints we are interested into *)
      let ints = A.init rand_bound (fun i -> i) in
      A.shuffle ~state:rng ints; (* random permute them *)
      ints
    ) seeds

(* compute k hashes *)
let hash pregen_rands idx2feat feat2acc_bound (dense_fp: dense): hashed =
  let k = A.length pregen_rands in
  let res = A.make k 0 in
  for i = 0 to k - 1 do
    let misses = ref 0 in
    let j = ref 0 in
    let rands = pregen_rands.(i) in
    let rand' = rands.(!j) in
    let test_feat_id = ref (idx2feat.(rand')) in
    let test_feat_val = ref (rand' - feat2acc_bound.(!test_feat_id)) in
    (* while is_red ... *)
    while !test_feat_val >= (BA1.unsafe_get dense_fp !test_feat_id) do
      incr misses; (* in the paper: Hashes[i]++ *)
      incr j;
      let rand = rands.(!j) in
      test_feat_id := idx2feat.(rand);
      test_feat_val := rand - feat2acc_bound.(!test_feat_id)
    done;
    A.unsafe_set res i !misses
  done;
  res

let estimate_jaccard (hash1: hashed) (hash2: hashed): float =
  let res = ref 0 in
  let k = A.length hash1 in
  for i = 0 to k - 1 do
    if (A.unsafe_get hash1 i) = (A.unsafe_get hash2 i) then
      incr res
  done;
  (float !res) /. (float k)

let estimate_distance (h1: hashed) (h2: hashed): float =
  1.0 -. (estimate_jaccard h1 h2)
