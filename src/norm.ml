
open Printf

module IntMap = BatMap.Int

(* max norm is probably to be preferred if we are going to minwise hash
 * the fingerprints later on *)
type norm = Max_norm (* max feature value in current instance *)
          | L1_norm (* Manhatan distance *)

let of_string = function
  | "l1" -> L1_norm
  | "max" -> Max_norm
  | other -> failwith (sprintf "Decoder: unknown norm: %s" other)

let map_norm style map =
  float
    (match style with
     | L1_norm -> IntMap.fold (fun _k v acc -> v + acc) map 0
     | Max_norm -> IntMap.fold (fun _k v acc -> max v acc) map 0)
