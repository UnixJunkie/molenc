(* counted Bloom filter *)

type t = int array array (* input feature index to output feature indexes mapping *)

(* n: input vector dimension
   k: number of "hash" functions;
      or number of output features turned on by a single input feature
   m: output vector dimension *)
let init n k m =
  let res = Array.make_matrix n k 0 in
  let rng = Random.State.make [|3141596|] in
  for i = 0 to n - 1 do
    for j = 0 to k - 1 do
      res.(i).(j) <- Random.State.int rng m
    done
  done;
  res
