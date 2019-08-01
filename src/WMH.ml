
(* Weighted Minwise Hashing in Constant Time

   Shrivastava, A. (2016).
   Simple and efficient weighted minwise hashing.
   In Advances in Neural Information Processing Systems (pp. 1498-1506). *)

module BA = Bigarray
module BA1 = BA.Array1
module Fp = Fingerprint

type dense = (int, BA.int8_unsigned_elt, BA.c_layout) BA1.t

type hashed = dense (* but the array will be shorter *)

let global_seed = 12345
let max_feat_val = 255
let k = 50

let rng = Random.State.make [|global_seed|]
let bound = (BatInt.pow 2 30) - 1
let seeds = Array.init k (fun _ -> BatRandom.State.int rng bound)

(* converse the sparse Fp.t type into a dense array of small positive ints *)
let to_dense max_feat fp =
  let res = BA1.create BA.int8_unsigned BA.C_layout max_feat in
  BA1.fill res 0;
  let n = BA1.dim fp in
  let i = ref 0 in
  while !i < n do
    (* unsafe *)
    let k1 = BA1.unsafe_get fp !i in
    let v1 = BA1.unsafe_get fp (!i + 1) in
    assert(k1 < max_feat && v1 < max_feature_value);
    let k1 = BA1.unsafe_get m1 !i in
  done;
  res
