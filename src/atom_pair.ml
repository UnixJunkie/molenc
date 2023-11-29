(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   An atom pair *)

(* FBR: TO DROP *)

type t = { src: PiEltHA.t; (* source atom *)
           dst: PiEltHA.t; (* destination atom *)
           dist: int } (* distance between them in bonds *)

(* canonicalization is obtained by sorting the types *)
let create x y dist =
  if PiEltHA.compare x y <= 0 then
    { src = x; dst = y; dist }
  else
    { src = y; dst = x; dist }

let to_string { src; dst; dist } =
  Printf.sprintf "%s-%d-%s" src dist dst

let dist x =
  x.dist
