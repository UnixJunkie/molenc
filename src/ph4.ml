
(* pharmacophore features supported by rdkit *)

open Printf

type t = Acc (* HB acceptor *)
       | Don (* HB donor *)
       | Pos (* pos. ionizable *)
       | Neg (* neg. ionizable *)
       | Hyd (* hydrohpobe *)
       | Lhy (* lumped hydrophobe *)
       | Znb (* Zn binder *)
       | Aro (* aromatic *)
       | Non (* none *)

let of_char = function
  | 'D' -> Don
  | 'A' -> Acc
  | 'P' -> Pos
  | 'N' -> Neg
  | 'a' -> Aro
  | 'H' -> Hyd
  | 'h' -> Lhy
  | 'Z' -> Znb
  | '_' -> Non
  | c -> failwith (sprintf "Ph4.of_char: unknown: %c" c)

let to_char = function
  | Don -> 'D'
  | Acc -> 'A'
  | Pos -> 'P'
  | Neg -> 'N'
  | Aro -> 'a'
  | Hyd -> 'H'
  | Lhy -> 'h'
  | Znb -> 'Z'
  | Non -> '_'

let to_string feat =
  String.make 1 (to_char feat)
