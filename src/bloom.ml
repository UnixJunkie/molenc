
(* A counted Bloom filter *)

module A = BatArray
module Log = Dolog.Log

type t = int array array (* input feature index (0..N-1) to output feature
                            indexes mapping (0..M-1) *)

(* n: input vector dimension
   k: number of "hash" functions;
      or number of output features "turned ON" by a single input feature
   m: output vector dimension *)
let init n k m =
  let res = Array.make_matrix n k 0 in
  let rng = Random.State.make [|3141596|] in
  for i = 0 to n - 1 do
    for j = 0 to k - 1 do
      res.(i).(j) <- Random.State.int rng m
    done
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
  res
