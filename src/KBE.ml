
(* KBE: "k-bits encoded" molecule *)

module A = Array
module Bst = Bst.Bisec_tree.Make (FpMol)
module Ht = Hashtbl
module L = List

type t = Bitv.t

let init mols =
  let k = L.length mols in
  let ht = Ht.create k in
  L.iteri (fun i mol ->
      (* vantage molecules must be unique *)
      assert(not (Ht.mem ht mol));
      Ht.add ht mol i
    ) mols;
  (* index them *)
  let bst = Bst.(create 1 Two_bands (A.of_list mols)) in
  (k, ht, bst)

let encode vmol2idx vpt mol =
  let k = Ht.length vmol2idx in
  let bits = Bitv.create k false in
  let nearby_mols = Bst.neighbors vpt 0.5 mol in
  L.iter (fun vmol ->
      let i = Ht.find vmol2idx vmol in
      Bitv.set bits i true
    ) nearby_mols;
  bits

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
