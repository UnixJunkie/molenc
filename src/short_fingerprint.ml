
(* Very short fingerprints, to encode an uncounted-unfolded FP into
   a bit vector (this gives a condensed/lower resolution version of the FP). *)

module BA = Bigarray
module BA1 = BA.Array1

type dense = WMH.dense
type t = Bitv.t

(* low resolution FP from an exact one *)
let of_dense nbits idx2feat feat2acc_bound (dense_fp: dense): t =
  let bits = Bitv.create nbits false in
  let step_size = (BA1.dim dense_fp) / nbits in
  let i = ref 0 in
  for j = 0 to nbits - 1 do
    let feat_id = idx2feat.(!i) in
    let feat_val = !i - feat2acc_bound.(feat_id) in
    if WMH.is_green dense_fp feat_id feat_val then
      Bitv.set bits j true;
    i := !i + step_size
  done;
  bits

(* |AnB| / |AuB| *)
let jaccard (hash1: t) (hash2: t): float =
  (float (Bitv.pop (Bitv.inter hash1 hash2))) /.
  (float (Bitv.pop (Bitv.union hash1 hash2)))

let distance (h1: t) (h2: t): float =
  1.0 -. (estimate_jaccard h1 h2)
