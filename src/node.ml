(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

module IntSet = BatSet.Int

type t = { typ: PiEltHA.t ; (* atom type *)
           succs: IntSet.t } (* indexes of its direct successors
                                in the molecular graph (it is bonded to them) *)

let create typ succs =
  { typ; succs }

let dummy = create PiEltHA.dummy IntSet.empty

let add_succ (n: t) (succ: int): t =
  create n.typ (IntSet.add succ n.succs)

let get_succs (n: t): IntSet.t =
  n.succs

let get_typ (n: t): PiEltHA.t =
  n.typ
