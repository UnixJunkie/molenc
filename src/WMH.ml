
(* Weighted Minwise Hashing in (amortized) Constant Time

   Shrivastava, A. (2016).
   Simple and efficient weighted minwise hashing.
   In Advances in Neural Information Processing Systems (pp. 1498-1506). *)

module A = Array
module BA = Bigarray
module BA1 = BA.Array1
module Fp = Fingerprint
module L = BatList

type dense = (int, BA.int8_unsigned_elt, BA.c_layout) BA1.t

(* type hashed = (int, BA.int16_unsigned_elt, BA.c_layout) BA1.t *)
type hashed = int array

(* let max_int16_unsigned = (BatInt.pow 2 16) - 1 *)

(* any feature value [x] must satisfy [x < feat_val_bound] *)
(* let feat_val_bound = 256 *)

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
    (* assert(v < feat_val_bound); *)
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
    (* assert(v < feat_val_bound); *)
    i := !i + 2
  done

(* compute the max value for each feature. *)
(* I.e. the columns' maximum if we put observations as rows
 * and features as columns in a data matrix *)
let bounds (max_feat_id: int) (train: Fp.t array): int array =
  let bounds = A.make max_feat_id 0 in
  A.iter (update_bounds bounds) train;
  bounds

(* in the paper, he defines is_green; but he samples until is_green becomes
 * true. It is more natural to sample while is_red *)
let is_red (arr: dense) (test_feat_id: int) (test_feat_val: int): bool =
  let feat_val = BA1.get arr test_feat_id in
  (feat_val = 0) || (test_feat_val > feat_val)

(* compute k hashes *)
let hash seeds bounds (dense_fp: dense): hashed =
  let k = A.length seeds in
  let feat_id_bound = BA1.dim dense_fp in
  (* let res = BA1.create BA.int16_unsigned BA.C_layout k in
   * BA1.fill res 0; *)
  let res = A.make k 0 in
  for i = 0 to k - 1 do
    let misses = ref 0 in
    let seed = A.get seeds i in
    let rng = Random.State.make [|seed|] in
    (* FBR: we could generate a single rand then modulo,
            if hashing speed really matters *)
    (* FBR: we could also generate enough rands in advance... *)
    let test_feat_id = ref (Random.State.int rng feat_id_bound) in
    let test_feat_val =
      let bound = 1 + bounds.(!test_feat_id) in
      ref (Random.State.int rng bound) in
    (* let test_feat_val = ref (Random.State.int rng feat_val_bound) in *)
    while is_red dense_fp !test_feat_id !test_feat_val do
      incr misses; (* Hashes[i]++ *)
      test_feat_id := Random.State.int rng feat_id_bound;
      test_feat_val :=
        let bound = 1 + bounds.(!test_feat_id) in
        Random.State.int rng bound
      (* test_feat_val := Random.State.int rng feat_val_bound *)
    done;
    (* assert(!misses <= max_int16_unsigned); (\* overflow? *\) *)
    (* BA1.set res i !misses *)
    res.(i) <- !misses
  done;
  res

let estimate_jaccard (hash1: hashed) (hash2: hashed): float =
  let res = ref 0 in
  (* let k = BA1.dim hash1 in
   * assert(BA1.dim hash2 = k); *)
  let k = A.length hash1 in
  assert(A.length hash2 = k);
  for i = 0 to k - 1 do
    (* if (BA1.get hash1 i) = (BA1.get hash2 i) then *)
    if hash1.(i) = hash2.(i) then
      incr res
  done;
  (float !res) /. (float k)

let estimate_distance (h1: hashed) (h2: hashed): float =
  1.0 -. (estimate_jaccard h1 h2)
