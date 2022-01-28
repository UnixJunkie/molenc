(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo, Japan.

   Counted atom pairs, with a new atom-typing scheme *)

open Printf

let () = Py.initialize ()

(* (atomic_num, nb_HA, nb_H, HA_used_val, formal_charge) cf. rdkit_wrapper.py *)
type atom = (int * int * int * int * int)
type pair = (atom * int * atom)

(* canonicalization is obtained by sorting the types *)
let create src dist dst =
  if compare src dst <= 0 then
    (src, dist, dst)
  else
    (dst, dist, src)

let string_of_atom (a,b,c,d,e) =
  sprintf "(%d,%d,%d,%d,%d)" a b c d e

let to_string (src, dist, dst) =
  Printf.sprintf "%s-%d-%s" (string_of_atom src) dist (string_of_atom dst)

let main () =
  failwith "not implemented yet"

let () = main ()
