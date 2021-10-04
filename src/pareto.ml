(* Copyright (C) 2021, Francois Berenger
   Tsuda laboratory,
   Graduate School of Frontier Sciences,
   The University of Tokyo,
   5-1-5 Kashiwa-no-ha,
   Kashiwa, Chiba 277-8561, Japan.

   Compute the Pareto front from a set of solutions.
   Higher scores (in each dimension) are considered better. *)

open Printf

module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

let parse_field_nums s =
  let strings = S.split_on_char ',' s in
  (* the user wants 1st field at index 1 (as in Awk), but on computers
     the 1st index is 0 *)
  L.map (fun num_str -> (int_of_string num_str) - 1) strings

type solution = float list

let sum_scores (s: solution): float =
  L.fsum s

let sum_scores_decr_sort solutions =
  L.sort (fun s1 s2 ->
      BatFloat.compare (sum_scores s2) (sum_scores s1)
    ) solutions

let is_dominated (s1: solution) (s2: solution): bool =
  L.for_all2 (fun p1 p2 -> p2 >= p1) s1 s2

let is_dominated_by_any candidate competitors =
  L.exists (is_dominated candidate) competitors

let remove_dominated candidate competitors =
  L.filter (fun competitor ->
      not (is_dominated competitor candidate)
    ) competitors

let pareto_front solutions =
  (* decr. sum scores sort *)
  let sorted = sum_scores_decr_sort solutions in
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
  let fields = S.split_on_char sep line in
  L.map (fun field_num ->
      let float_str = L.at fields field_num in
      try (* robust float parsing *)
        Scanf.sscanf float_str "%f" (fun x -> x)
      with exn ->
        let () =
          Log.fatal "Pareto.solution_of_string: cannot parse: '%s' in '%s'"
            float_str line in
        raise exn
    ) field_nums

let print_solution out sol =
  L.iteri (fun i x ->
      fprintf out (if i > 0 then " %f" else "%f") x
    ) sol;
  fprintf out "\n"

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file\n  \
               -d <char>: field separator (default=\\t)\n  \
               -f <int>,<int>[,...]: fields to filter on \
               (first field index=1; same as awk)\n" Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let all_lines_uniq =
    L.unique_cmp ~cmp:S.compare
      (* ignore comment lines / starting with '#' *)
      (LO.filter input_fn (fun l -> not (S.starts_with l "#"))) in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field_nums_str = CLI.get_string ["-f"] args in
  CLI.finalize(); (* ------------------------------------------------------- *)
  let field_nums = parse_field_nums field_nums_str in
  let solutions = L.map (solution_of_string sep field_nums) all_lines_uniq in
  let front = pareto_front solutions in
  L.iter (print_solution stdout) front

let () = main ()
