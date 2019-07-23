(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

open Printf

type t = Single of int (* encoding radius *)
       | Multi of int * int (* start-stop encoding radii *)

let of_string s =
  if BatString.contains s ':' then
    let istr, jstr = BatString.split s ~by:":" in
    let i, j = int_of_string istr, int_of_string jstr in
    assert(i >= 0 && j >= 0 && i <= j);
    Multi (i, j)
  else
    Single (int_of_string s)

let to_string = function
  | Single i -> sprintf "%d" i
  | Multi (i, j) -> sprintf "%d..%d" i j

let to_list = function
  | Single i -> [i]
  | Multi (i, j) -> BatList.range i `To j
