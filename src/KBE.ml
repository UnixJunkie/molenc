
(* KBE: "k-bits encoded" molecule *)

module A = Array
module Bst = Bst.Bisec_tree.Make (FpMol)
module Ht = Hashtbl
module L = List

type t = Bitv.t

let init mols =
  let k = L.length mols in
  let ht = Ht.create k in
  let thresholds =
    let rng = Random.State.make_self_init () in
    A.init k (fun _ -> Random.State.float rng 1.0) in
  L.iteri (fun i mol ->
      (* vantage molecules must be unique *)
      assert(not (Ht.mem ht mol));
      Ht.add ht mol i
    ) mols;
  (* index them *)
  let bst = Bst.(create 1 Two_bands (A.of_list mols)) in
  (ht, thresholds, bst)

let encode vmol2idx thresholds bst mol =
  let k = Ht.length vmol2idx in
  let bits = Bitv.create k false in
  for i = 0 to k - 1 do
    let t = thresholds.(i) in
    let nearby_mols = Bst.neighbors mol t bst in
    L.iter (fun vmol ->
        let j = Ht.find vmol2idx vmol in
        if j = i then
          Bitv.set bits i true
      ) nearby_mols
  done;
  (* let nearby_mols = Bst.neighbors mol 0.5 bst in
   * L.iter (fun vmol ->
   *     let i = Ht.find vmol2idx vmol in
   *     Bitv.set bits i true
   *   ) nearby_mols; *)
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
