
(* small module to help compute atom environments from *.ph4 files *)

open Printf

type t = { idx: int;
           typ: Ph4.t }

let create (idx: int) (typ: Ph4.t): t =
  (* indexes start at 1 in MOL2 files *)
  { idx = idx - 1; typ }

let dummy = create (-1) Ph4.Non

let of_ph4_line (l: string): t =
  try Scanf.sscanf l "%d %c"
        (fun idx char -> create idx (Ph4.of_char char))
  with Scanf.Scan_failure msg ->
    failwith (sprintf "Ph4_atom.of_ph4_line: could not parse: %s: %s" l msg)

let to_string (a: t): string =
  sprintf "%d %s" a.idx (Ph4.to_string a.typ)
