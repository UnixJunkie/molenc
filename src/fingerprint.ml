(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

module A = Array
module BA = Bigarray
module BA1 = BA.Array1
module Ht = Hashtbl
module IntMap = BatMap.Int
module L = MyList

(* an unfolded-counted FP *)
(* the int16 bigarray trick reduces memory consumption by four times compared
 * to regular 64 bits OCaml integers; thanks to Oleg for suggesting it
 * and to Chet Murthy for suggesting arrays *)
type t = (int, BA.int_elt, BA.c_layout) BA1.t

let create_BA1 n =
  BA1.create BA.Int BA.C_layout n

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
             Utls.enforce_f (k >= 0 && k > !previous && v > 0)
               (fun () -> "Fingerprint.of_string: invalid line: " ^ s);
             previous := k;
             incr n;
             (k, v)
          )
      ) s in
  let res = create_BA1 (2 * !n) in
  let i = ref 0 in
  L.iter (fun (k, v) ->
      BA1.unsafe_set res !i k;
      incr i;
      BA1.unsafe_set res !i v;
      incr i
    ) kvs;
  res

let to_string (x: t): string =
  let buff = Buffer.create 80 in
  let n = BA1.dim x in
  let i = ref 0 in
  while !i < n do
    let k = BA1.unsafe_get x !i in
    let v = BA1.unsafe_get x (!i + 1) in
    if !i = 0 then
      Printf.bprintf buff "%d:%d" k v
    else
      Printf.bprintf buff ";%d:%d" k v;
    i := !i + 2
  done;
  Buffer.contents buff

let max_feat_id x =
  let n = BA1.dim x in
  BA1.get x (n - 2)

let nb_features x =
  1 + (max_feat_id x)

(* sparse to dense conversion *)
let to_dense (max_len: int) (x: t): int array =
  let res = A.make max_len 0 in
  let n = BA1.dim x in
  let i = ref 0 in
  while !i < n do
    res.(BA1.unsafe_get x !i) <- BA1.unsafe_get x (!i + 1);
    i := !i + 2
  done;
  res

(* sparse to dense printf;
   without creation of the intermediate dense array (for perfs) *)
let to_dense_printf nb_features (x: t): unit =
  let n = BA1.dim x in
  let i = ref 0 in
  let j = ref 0 in
  while !i < n do
    let k = BA1.unsafe_get x !i in
    while !j < k do
      Printf.printf " 0";
      incr j
    done;
    let v = BA1.unsafe_get x (!i + 1) in
    Printf.printf " %d" v;
    incr j;
    i := !i + 2
  done;
  while !j < nb_features do
    Printf.printf " 0";
    incr j
  done

let sum_min_max (m1: t) (m2: t): (int * int) =
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
  (!icard, !ucard)

(* tani(A,B) = |inter(A,B)| / |union(A,B)|
             = sum(min_i) / sum(max_i) *)
let tanimoto (m1: t) (m2: t): float =
  let icard, ucard = sum_min_max m1 m2 in
  if ucard = 0 then 0.0
  else (float icard) /. (float ucard)

(* tanimoto distance (this _is_ a metric) *)
let distance x y =
  1.0 -. (tanimoto x y)

(* convert to int map: feat_id -> feat_val *)
let key_values fp =
  let res = ref IntMap.empty in
  let len = BA1.dim fp in
  let i = ref 0 in
  while !i < len do
    let k = BA1.unsafe_get fp !i in
    let v = BA1.unsafe_get fp (!i + 1) in
    res := IntMap.add k v !res;
    i := !i + 2
  done;
  !res

let key_value_pairs fp =
  let res = ref [] in
  let len = BA1.dim fp in
  let i = ref 0 in
  while !i < len do
    let k = BA1.unsafe_get fp !i in
    let v = BA1.unsafe_get fp (!i + 1) in
    res := (k, v) :: !res;
    i := !i + 2
  done;
  !res

(* iterate given function on all key-value pairs *)
let kv_iter f fp =
  let len = BA1.dim fp in
  let i = ref 0 in
  while !i < len do
    f (BA1.unsafe_get fp !i) (BA1.unsafe_get fp (!i + 1));
    i := !i + 2
  done

let drop_features to_drop fp =
  let kept =
    let kvs = key_values fp in
    IntMap.filter (fun k _v ->
        not (Ht.mem to_drop k)
      ) kvs in
  let n = IntMap.cardinal kept in
  let res = create_BA1 (2 * n) in
  let i = ref 0 in
  IntMap.iter (fun k v ->
      BA1.unsafe_set res !i k;
      incr i;
      BA1.unsafe_set res !i v;
      incr i
    ) kept;
  res

let sum_values fp =
  let len = BA1.dim fp in
  let i = ref 1 in (* values start at 1 *)
  let total = ref 0 in
  while !i < len do
    total := !total + (BA1.unsafe_get fp !i);
    i := !i + 2
  done;
  !total
