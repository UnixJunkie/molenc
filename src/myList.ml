(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

include BatList

let to_string to_str l =
  let buff = Buffer.create 80 in
  Buffer.add_char buff '[';
  iteri (fun i x ->
      if i > 0 then Buffer.add_char buff ';';
      Buffer.add_string buff (to_str x);
    ) l;
  Buffer.add_char buff ']';
  Buffer.contents buff

let of_string of_str s =
  let s' = BatString.chop ~l:1 ~r:1 s in
  if s' = "" then
    (* the empty list case was not handled in the past *)
    []
  else
    begin
      if String.contains s' ']' then
        failwith ("MyList.of_string: sub lists inside: " ^ s);
      map of_str (BatString.split_on_string s' ~by:";")
    end

(* only map 'f' on elements satisfying 'p' *)
let filter_map p f l =
  let res =
    fold_left (fun acc x ->
        if p x then (f x) :: acc
        else acc
      ) [] l in
  rev res

(* split a list into n parts (the last part might have
   a different number of elements) *)
let nparts n l =
  let len = length l in
  let res = ref [] in
  let curr = ref l in
  let m = BatFloat.round_to_int (float len /. float n) in
  for _ = 1 to n - 1 do
    let xs, ys = takedrop m !curr in
    curr := ys;
    res := xs :: !res
  done;
  rev (!curr :: !res)

(* create folds of cross validation; each fold consists in (train, test) *)
let cv_folds n l =
  let test_sets = nparts n l in
  let rec loop acc prev curr =
    match curr with
    | [] -> acc
    | x :: xs ->
      let before_after = flatten (rev_append prev xs) in
      let prev' = x :: prev in
      let train_test = (before_after, x) in
      let acc' = train_test :: acc in
      loop acc' prev' xs in
  loop [] [] test_sets

(* List.combine for 4 lists *)
let combine4 l1 l2 l3 l4 =
  let rec loop acc = function
    | ([], [], [], []) -> rev acc
    | (w :: ws, x :: xs, y :: ys, z :: zs) ->
      loop ((w, x, y, z) :: acc) (ws, xs, ys, zs)
    | _ -> raise (Invalid_argument "MyList.combine4: list lengths differ")
  in
  loop [] (l1, l2, l3, l4)

let really_take n l =
  let res = take n l in
  assert(length res = n);
  res

(* non reproducible randomization of a list *)
let random_shuffle l =
  let rng = BatRandom.State.make_self_init () in
  shuffle ~state:rng l

let rev_combine l1 l2 =
  let rec loop acc l r =
    match (l, r) with
    | ([], []) -> acc
    | (x :: xs, y :: ys) -> loop ((x, y) :: acc) xs ys
    | _ -> raise (Invalid_argument "MyList.rev_combine: list lengths differ")
  in
  loop [] l1 l2

(* filter using bit-mask [m] *)
let filter_mask m l =
  let rec loop acc = function
    | [] -> acc
    | (p, x) :: rest -> loop (if p then x :: acc else acc) rest
  in
  loop [] (rev_combine m l)
