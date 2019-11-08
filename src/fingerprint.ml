(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

module A = Array
module BA = Bigarray
module BA1 = BA.Array1
module IntMap = BatMap.Int
module L = MyList

(* an unfolded-counted FP *)
(* the int16 bigarray trick reduces memory consumption by four times compared
 * to regular 64 bits OCaml integers; thanks to Oleg for suggesting it
 * and to Chet Murthy for suggesting arrays *)
type t = (int, BA.int_elt, BA.c_layout) BA1.t

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
  let res = BA1.create BA.Int BA.C_layout (2 * !n) in
  let i = ref 0 in
  L.iter (fun (k, v) ->
      BA1.unsafe_set res !i k;
      incr i;
      BA1.unsafe_set res !i v;
      incr i
    ) kvs;
  res

let nb_features x =
  let n = BA1.dim x in
  1 + (BA1.get x (n - 2))

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

(* (\* Euclidian distance; in case we work with 2D points instead of molecules
 *    HACK: POINTS HAVE COORDINATES AS INT BUT THEY NEED TO BE DIVIDED BY 1000
 *          SINCE THEY ARE SUPPOSED TO BE 2D POINTS IN THE [0,1],[0,1] QUADRANT *\)
 * let distance_2D m1 m2 =
 *   let len1 = BA1.dim m1 in
 *   let len2 = BA1.dim m2 in
 *   assert(len1 = 4 && len2 = 4);
 *   let ix1 = BA1.unsafe_get m1 0 in
 *   let x1  = (float (BA1.unsafe_get m1 1)) /. 1000.0 in
 *   let iy1 = BA1.unsafe_get m1 2 in
 *   let y1  = (float (BA1.unsafe_get m1 3)) /. 1000.0 in
 *   assert(ix1 = 0 && iy1 = 1);
 *   let ix2 = BA1.unsafe_get m2 0 in
 *   let x2  = (float (BA1.unsafe_get m2 1)) /. 1000.0 in
 *   let iy2 = BA1.unsafe_get m2 2 in
 *   let y2  = (float (BA1.unsafe_get m2 3)) /. 1000.0 in
 *   assert(ix2 = 0 && iy2 = 1);
 *   let dx = x1 -. x2 in
 *   let dy = y1 -. y2 in
 *   sqrt (dx *. dx +. dy *. dy) *)
