(* Copyright (C) 2021, Francois Berenger
   Tsuda laboratory,
   Graduate School of Frontier Sciences,
   The University of Tokyo,
   5-1-5 Kashiwa-no-ha,
   Kashiwa, Chiba 277-8561, Japan.

   Compute the Pareto front from a set of solutions.
   Higher scores (in each dimension) are considered better. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

let parse_field_nums s =
  let strings = S.split_on_char ',' s in
  (* the user wants 1st field at index 1 (as in Awk), but on computers
     the 1st index is 0 *)
  L.map (fun num_str -> (int_of_string num_str) - 1) strings

type solution = float array

type ranked_solution = int array

let sum_ranks (s: ranked_solution): int =
  A.sum s

(* array of ranks for given float array's sorted copy *)
let rank_array (a: float array): (float array * int array) =
  let sorted = A.copy a in
  (* sort by decreasing scores: highest score gets lowest rank *)
  A.sort (fun x y -> BatFloat.compare y x) sorted;
  let ranks =
    let n = A.length a in
    A.make n 0 in
  let curr_rank = ref (-1) in
  let previous = ref neg_infinity in
  A.iteri (fun i x ->
      if x <> !previous then
        (incr curr_rank;
         previous := x;
         ranks.(i) <- !curr_rank)
      else
        ranks.(i) <- !curr_rank
    ) sorted;
  (sorted, ranks)

(* unit test for rank_array *)
let () =
  assert(snd (rank_array [|-1.; 0.; 0.; 1.; 1.; 2.; 2.|]) =
         [|0; 0; 1; 1; 2; 2; 3|])

let score2rank_ht (scores: float array): (float, int) Ht.t =
  let sorted_scores, corresp_ranks = rank_array scores in
  let ht = Ht.create (A.length sorted_scores) in
  A.iter2 (fun score rank ->
      if not (Ht.mem ht score) then
        Ht.add ht score rank
    ) sorted_scores corresp_ranks;
  ht

let extract_coord (i: int) (sols: solution list): float array =
  let n = L.length sols in
  let res = A.make n 0.0 in
  L.iteri (fun j sol ->
      res.(j) <- sol.(i)
    ) sols;
  res

let rank_solutions (sols: solution list): ranked_solution list =
  let ndims = A.length (L.hd sols) in
  (* create all score2rank hts *)
  let score2rank_hts =
    A.init ndims (fun coord ->
        let xs = extract_coord coord sols in
        score2rank_ht xs
      ) in
  (* apply them on each coordinate *)
  L.map (fun sol ->
      A.mapi (fun i x ->
          let ht = score2rank_hts.(i) in
          Ht.find ht x
        ) sol
    ) sols

let ranked_sol2sol_ht sols ranked_sols =
  let ht = Ht.create (L.length sols) in
  L.iter2 (fun sol ranked_sol ->
      Ht.replace ht ranked_sol sol
    ) sols ranked_sols;
  ht

(* rank by increasing sum of ranks *)
let sum_ranks_sort (sols: ranked_solution list): ranked_solution list =
  L.sort (fun s1 s2 ->
      BatInt.compare (sum_ranks s1) (sum_ranks s2)
    ) sols

let is_dominated (s1: ranked_solution) (s2: ranked_solution): bool =
  A.for_all2 (fun r1 r2 -> r2 <= r1) s1 s2

let is_dominated_by_any candidate competitors =
  L.exists (is_dominated candidate) competitors

let remove_dominated candidate competitors =
  L.filter (fun competitor ->
      not (is_dominated competitor candidate)
    ) competitors

let pareto_front solutions =
  (* Heuristic: try to put good solutions at the front, so they
     will prune weaker solutions and hence lower algorithmic complexity
     (a good solution has a low sum of ranks) *)
  let sorted = sum_ranks_sort solutions in
  (* loop until no dominated solution is left *)
  let rec loop acc = function
    | [] -> L.rev acc
    | x :: xs ->
      let xs' = remove_dominated x xs in
      let acc' =
        if is_dominated_by_any x xs' then
          acc
        else
          x :: acc in
      loop acc' xs' in
  loop [] sorted

let solution_of_string sep field_nums line =
  let num_fields = L.length field_nums in
  let res = A.make num_fields 0.0 in
  let fields = S.split_on_char sep line in
  L.iteri (fun i field_num ->
      let float_str = L.at fields field_num in
      let x =
        try (* robust float parsing *)
          Scanf.sscanf float_str "%f" (fun x -> x)
        with exn ->
          let () =
            Log.fatal "Pareto.solution_of_string: cannot parse: '%s' in '%s'"
              float_str line in
          raise exn in
      res.(i) <- x
    ) field_nums;
  res

let string_of_solution sol =
  let buff = Buffer.create 80 in
  A.iteri (fun i x ->
      if i > 0 then Buffer.add_char buff ' ';
      Buffer.add_string buff (string_of_float x)
    ) sol;
  Buffer.contents buff

let maybe_create_preserve_ht sep field_nums preserve_line all_lines =
  if not preserve_line then
    None
  else
    let ht = Ht.create (L.length all_lines) in
    L.iter (fun line ->
        let sol = solution_of_string sep field_nums line in
        Ht.replace ht sol line
      ) all_lines;
    Some ht

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename>: input file\n  \
              -d <char>: field separator (default=\t)\n  \
              -f <int>,<int>[,...]: fields to filter on \
              (first field index=1; same as awk)\n  \
              [--preserve]: preserve whole input line in output\n"
       Sys.argv.(0);
     exit 1
    );
  let input_fn = CLI.get_string ["-i"] args in
  let all_lines =
    (* ignore comments / lines starting with '#' *)
    LO.filter input_fn (fun l -> not (S.starts_with l "#")) in
  let all_lines_uniq = L.unique_cmp ~cmp:S.compare all_lines in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field_nums_str = CLI.get_string ["-f"] args in
  let preserve_line = CLI.get_set_bool ["--preserve"] args in
  CLI.finalize(); (* ------------------------------------------------------- *)
  let field_nums = parse_field_nums field_nums_str in
  let maybe_sol2preserved =
    maybe_create_preserve_ht sep field_nums preserve_line all_lines in
  let solutions = L.map (solution_of_string sep field_nums) all_lines_uniq in
  let ranked_solutions = rank_solutions solutions in
  let ranked_sol2sol = ranked_sol2sol_ht solutions ranked_solutions in
  let front = pareto_front ranked_solutions in
  match maybe_sol2preserved with
  | None ->
    L.iter (fun ranked_sol ->
        let sol = Ht.find ranked_sol2sol ranked_sol in
        fprintf stdout "%s\n" (string_of_solution sol)
      ) front
  | Some ht ->
    (* Pareto front coordinates then preserved line *)
    L.iter (fun ranked_sol ->
        let sol = Ht.find ranked_sol2sol ranked_sol in
        let preserved = Ht.find ht sol in
        fprintf stdout "%s %s\n" (string_of_solution sol) preserved
      ) front

let () = main ()
