(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* very small atom module to help compute atom environments (a la molprint2d)
   from MOL2 files *)

open Printf

type t = { idx: int;
           typ: Sybyl.t }

let create (idx: int) (typ: Sybyl.t): t =
  (* indexes start at 1 in MOL2 files *)
  { idx = idx - 1; typ }

let dummy = create (-1) Sybyl.Du

(* example line (output of OpenBabel):
   "      1 S          -0.0218    1.7554    0.0117 S.2     1  LIG1       -0.0637" *)
let of_mol2_line (l: string): t =
  try Scanf.sscanf l " %d %s %f %f %f %s@ %s@ %s@ %s"
        (fun idx _name _x _y _z typ _bs0 _bs1 _bs2 ->
           create idx (Sybyl.of_string typ))
  with Scanf.Scan_failure msg ->
    failwith (sprintf "Atom.of_mol2_line: could not parse: %s: %s" l msg)

let to_string (a: t): string =
  sprintf "%d %s" a.idx (Sybyl.to_string a.typ)
