(* Copyright (C) 2020, Francois Berenger

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

(* Example first line: ^#radius=0..1$
   i.e. the radius indicator is everything after "#radius=" *)
let of_dictionary_header fn =
  let header = Utls.get_first_line fn in
  let prfx = "#radius=" in
  let prfx_len = String.length prfx in
  Utls.enforce (BatString.starts_with header prfx)
    "Scale.of_dictionary_header: not a circular FP dictionary header; \
     --pairs CLI option probably missing";
  let s = BatString.lchop ~n:prfx_len header in
  of_string s

let to_string = function
  | Single i -> sprintf "%d" i
  | Multi (i, j) -> sprintf "%d:%d" i j

let to_list = function
  | Single i -> [i]
  | Multi (i, j) -> BatList.range i `To j
