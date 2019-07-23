(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* very small bond module to help compute atom environments (a la molprint2d)
   from MOL2 files *)

open Printf

type t = { idx: int ;
           src: int ;
           dst: int }

(* indexes start at 1 in the MOL2 file format *)
let create (idx: int) (src: int) (dst: int): t =
  { idx = idx - 1; src = src - 1; dst = dst - 1 }

let dummy = create (-1) (-1) (-1)

(* example line (output of OpenBabel):
   "     1     1    11   ar" *)
let of_mol2_line l =
  try Scanf.sscanf l " %d %d %d %s"
        (fun idx src dst _bullshit -> create idx src dst)
  with _ -> failwith ("Bond.of_mol2_line: could not parse: " ^ l)

let to_string (a: t): string =
  sprintf "%d %d %d" a.idx a.src a.dst
