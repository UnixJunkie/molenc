
(* Very short fingerprints, to encode an uncounted-unfolded FP into
   a bit vector (this gives a condensed/lower resolution version of the FP). *)

module A = Array
module BA = Bigarray
module BA1 = BA.Array1

type dense = WMH.dense
type t = Bitv.t

(* low resolution FP from an exact one *)
let of_dense nbits idx2feat feat2acc_bound (dense_fp: dense): t =
  let bits = Bitv.create nbits false in
  (* Log.debug "total length: %d" (A.length idx2feat);
   * Log.debug "nbits: %d" nbits; *)
  let step_size = (A.length idx2feat) / nbits in
  (* Log.debug "step size: %d" step_size; *)
  let i = ref 0 in
  for j = 0 to nbits - 1 do
    let feat_id = idx2feat.(!i) in
    let feat_val = !i - feat2acc_bound.(feat_id) in
    if WMH.is_green dense_fp feat_id feat_val then
      Bitv.set bits j true;
    i := !i + step_size
  done;
  bits

(* string of 0s and 1s *)
let to_string x =
  let res = Bytes.create (Bitv.length x) in
  Bitv.iteri (fun i b ->
      Bytes.set res i (if b then '1' else '0')
    ) x;
  Bytes.to_string res

(* |AnB| / |AuB| *)
let jaccard (hash1: t) (hash2: t): float =
  (float (Bitv.pop (Bitv.bw_and hash1 hash2))) /.
  (float (Bitv.pop (Bitv.bw_or hash1 hash2)))

let distance (h1: t) (h2: t): float =
  1.0 -. (jaccard h1 h2)
