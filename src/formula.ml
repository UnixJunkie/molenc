
module A = BatArray
module L = BatList
module Log = Dolog.Log

type formula_item = Element of string
                  | Count of int

(* parse list of tokens *)
let rec count_elements = function
  | [] -> []
  | [Element symb] -> [(symb, 1)]
  | [Count _] -> assert(false) (* should have been processed before *)
  | (Element s0) :: (Element s1) :: rest ->
    (s0, 1) :: count_elements (Element s1 :: rest)
  | (Element symb) :: (Count c) :: rest ->
    (symb, c) :: (count_elements rest)
  | _ -> assert(false)

let parse_int s =
  try int_of_string s
  with exn ->
    (Log.fatal "Formula.parse_int: cannot parse: %s" s;
     raise exn)

(* formula -> int *)
let encode _debug f =
  let element_counts = A.make 119 0 in
  (* lexer: tokenize chemical elements starting from two chars ones *)
  let element_counts_0 =
    Str.bounded_full_split Ptable.elements_regexp f 1024 in
  let element_counts_1 =
    L.map (function Str.Delim symbol -> Element symbol
                  | Str.Text count -> Count (parse_int count)
      ) element_counts_0 in
  let element_counts_2 = count_elements element_counts_1 in
  L.iter (fun (symb, count) ->
      let anum = Ptable.anum_of_symbol symb in
      (* robust even to extended formulas like CH3CH2OH
         instead of the proper C2H6O *)
      element_counts.(anum) <- element_counts.(anum) + count
    ) element_counts_2;
  (* potentially too large number to fit OCaml's 64 bits signed integers *)
  let big_int =
    A.fold_lefti (fun acc anum count ->
        if count > 0 then
          let p = Z.of_int (Ptable.prime_for_anum anum) in
          (* sum of powers of primes; also called "Godel numbering" *)
          Z.mul acc (Z.pow p count)
        else
          acc
      ) Z.one element_counts in
  Z.to_int big_int

let compute_exponent composite prime =
  let rec loop acc x =
    if x mod prime = 0 then
      loop (acc + 1) (x / prime)
    else
      acc in
  loop 0 composite

(* int -> formula *)
let decode _debug _x =
  (* prime x is a factor if y mod x = 0 *)
  failwith "FBR: TO IMPLEMENT"
