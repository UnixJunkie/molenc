
(* KBE: "k-bits encoded" molecule *)

module A = Array
module Bst = Bst.Bisec_tree.Make (FpMol)
module Ht = Hashtbl
module L = List

type t = Bitv.t

let init rng mols =
  let k = L.length mols in
  let ht = Ht.create k in
  let thresholds = A.init k (fun _ -> Random.State.float rng 1.0) in
  (* to guarantee each bit adds some information:
     vantage molecules must be unique *)
  L.iter (fun mol ->
      let fp = FpMol.get_fp mol in
      assert(not (Ht.mem ht fp));
      Ht.add ht fp ()
    ) mols;
  (A.of_list mols, thresholds)

let encode vmols thresholds mol =
  let k = A.length vmols in
  let bits = Bitv.create k false in
  for i = 0 to k - 1 do
    if FpMol.dist mol (A.unsafe_get vmols i) <= (A.unsafe_get thresholds i) then
      Bitv.set bits i true
  done;
  (FpMol.get_name mol, bits)

let to_string x =
  Bitv.L.to_string x

(* This is the exact Tanimoto/Jaccard formula,
 * but this FP is a kind of hash of the initial fingerprint,
 * so let's call this Tanimoto an estimate *)
let estimate_tanimoto x y =
  (float (Bitv.pop (Bitv.bw_and x y))) /.
  (float (Bitv.pop (Bitv.bw_or x y)))

let estimate_distance x y =
  1.0 -. (estimate_tanimoto x y)
