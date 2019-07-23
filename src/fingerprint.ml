(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* open Printf *)

module A = Array
module BA = Bigarray
module BA1 = BA.Array1
module IntMap = BatMap.Int
module L = MyList

(* an unfolded-counted FP *)
(* the int16 bigarray trick reduces memory consumption by four times compared
 * to regular 64 bits OCaml integers; thanks to Oleg for suggesting it
 * and to Chet Murthy for suggesting arrays *)
type t = (int, BA.int16_unsigned_elt, BA.c_layout) BA1.t

let of_string s: t =
  let previous = ref (-1) in
  let n = ref 0 in
  let kvs =
    L.of_string (fun str ->
        Scanf.sscanf str "%d:%d"
          (fun k v ->
             (* indices are >= 0 *)
             (* indices are incr. sorted *)
             (* feature counts are > 0 *)
             (* keys and values all fit in an unsigned int 16 *)
             assert(k >= 0 && k > !previous && v > 0 &&
                    k <= 65535 && v <= 65535);
             previous := k;
             incr n;
             (k, v)
          )
      ) s in
  let res = BA1.create BA.Int16_unsigned BA.C_layout (2 * !n) in
  let i = ref 0 in
  L.iter (fun (k, v) ->
      BA1.unsafe_set res !i k;
      incr i;
      BA1.unsafe_set res !i v;
      incr i
    ) kvs;
  res

(* tani(A,B) = |inter(A,B)| / |union(A,B)|
             = sum(min_i) / sum(max_i) *)
let tanimoto (m1: t) (m2: t): float =
  let icard = ref 0 in
  let ucard = ref 0 in
  let len1 = BA1.dim m1 in
  let len2 = BA1.dim m2 in
  let i = ref 0 in
  let j = ref 0 in
  while !i < len1 && !j < len2 do
    (* unsafe *)
    let k1 = BA1.unsafe_get m1 !i in
    let v1 = BA1.unsafe_get m1 (!i + 1) in
    let k2 = BA1.unsafe_get m2 !j in
    let v2 = BA1.unsafe_get m2 (!j + 1) in
    (* process keys in increasing order *)
    if k1 < k2 then
      (ucard := !ucard + v1;
       i := !i + 2)
    else if k2 < k1 then
      (ucard := !ucard + v2;
       j := !j + 2)
    else (* k1 = k2 *)
    if v1 <= v2 then
      (icard := !icard + v1;
       ucard := !ucard + v2;
       i := !i + 2;
       j := !j + 2)
    else
      (icard := !icard + v2;
       ucard := !ucard + v1;
       i := !i + 2;
       j := !j + 2)
  done;
  incr i; (* go to value *)
  while !i < len1 do (* finish m1; unsafe *)
    ucard := !ucard + (BA1.unsafe_get m1 !i);
    i := !i + 2
  done;
  incr j; (* go to value *)
  while !j < len2 do (* finish m2; unsafe *)
    ucard := !ucard + (BA1.unsafe_get m2 !j);
    j := !j + 2
  done;
  if !ucard = 0 then 0.0
  else (float !icard) /. (float !ucard)

(* tanimoto distance (this _is_ a metric) *)
let distance x y =
  1.0 -. (tanimoto x y)
