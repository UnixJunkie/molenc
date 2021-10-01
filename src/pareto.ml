
(* Compute the Pareto front for a set of solutions *)

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
  let _all_lines_uniq =
    L.unique_cmp ~cmp:S.compare
      (* ignore comment lines (starting with '#') *)
      (LO.filter input_fn (fun l -> not (S.starts_with l "#"))) in
  let _sep = CLI.get_char_def ["-d"] args '\t' in
  let field_nums_str = CLI.get_string ["-f"] args in
  CLI.finalize();
  let _field_nums = parse_field_nums field_nums_str in
  failwith "not implemented yet"

let () = main ()
