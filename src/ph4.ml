
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

let of_char = function
  | 'D' -> Don
  | 'A' -> Acc
  | 'P' -> Pos
  | 'N' -> Neg
  | 'a' -> Aro
  | 'H' -> Hyd
  | 'h' -> Lhy
  | 'Z' -> Znb
  | c -> failwith (sprintf "Ph4.of_char: unknown: %c" c)
