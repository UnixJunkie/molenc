
(* A counted Bloom filter *)

module A = BatArray
module Fp = Fingerprint
module L = BatList
module Log = Dolog.Log

type t = int array array (* input feature index (0..N-1) to output feature
                            indexes mapping (0..M-1) *)

let distinct_rands rng n bound =
  let rec loop acc count =
    if count = n then
      acc
    else
      let cand = Random.State.int rng bound in
      if List.mem cand acc then
        loop acc count (* retry *)
      else
        loop (cand :: acc) (count + 1) in
  loop [] 0

(* n: input vector dimension
   k: number of "hash" functions;
      number of output features "turned ON" by a single input feature
   m: output vector dimension *)
let init n k m =
  let res = Array.make_matrix n k 0 in
  let rng = Random.State.make [|3141596|] in
  for i = 0 to n - 1 do
    let rands = distinct_rands rng k m in
    L.iteri (fun j rand ->
        res.(i).(j) <- rand
      ) rands
  done;
  (* log the number of collisions
     (different input features mapping to the same set of output features *)
  let collisions = ref 0 in
  let sorted = A.copy res in
  A.sort compare sorted;
  for i = 1 to n - 1 do
    if sorted.(i - 1) = sorted.(i) then
      incr collisions;
  done;
  (if !collisions > 0 then
     Log.warn "Bloom.init(%d,%d,%d): %d collisions" n k m !collisions
  );
  (n, k, m, res)

let encode (_n, k, m, mappings) fp =
  let kvs = Fp.key_value_pairs fp in (* sparse input vector *)
  let res = A.make m 0 in (* dense output vector *)
  L.iter (fun (key, value) ->
      let output_indexes = mappings.(key) in
      (* increment all corresponding output features *)
      for i = 0 to k - 1 do
        let j = output_indexes.(i) in
        res.(j) <- res.(j) + value
      done
    ) kvs;
  res
