
(* merge score files using average of normalized scores
   (an unweighted score merging method) *)

open Printf

module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module String = BatString
module StringSet = BatSet.String
module Utls = Molenc.Utls

type score_spec = { fn: string;
                    sep: char;
                    name_field: int;
                    score_field: int;
                    increasing: bool }

let parse_spec sep str =
  try
    Scanf.sscanf str "%s@:%d:%d" (fun fn name score ->
        { fn;
          sep;
          name_field = name - 1;
          score_field = (abs score) - 1;
          increasing = score < 0 }
      )
  with exn -> (Log.fatal "Merge.parse_spec: cannot parse: %s" str;
               raise exn)

let parse_ifs sep s =
  let specs =
    try String.split_on_string s ~by:","
    with Not_found -> [s] in
  L.map (parse_spec sep) specs

let normalize_scores spec =
  let name_raw_scores =
    Utls.maybe_map_on_lines_of_file spec.fn (fun line ->
        try
          let name = String.cut_on_char spec.sep spec.name_field line in
          let score_field = String.cut_on_char spec.sep spec.score_field line in
          let score = Scanf.sscanf score_field "%f" (fun x -> x) in
          Some (name, score)
        with _exn -> None
      ) in
  let raw_scores = L.map snd name_raw_scores in
  let mini, maxi = L.min_max raw_scores in
  let avg = Utls.faverage raw_scores in
  let std = Utls.stddev raw_scores in
  Log.info "fn: %s%s min:avg+/-std:max: %g:%g+/-%g:%g"
    spec.fn (if spec.increasing then "(incr)" else "")
    mini avg std maxi;
  let sign = if spec.increasing then -1.0 else +1.0 in
  L.map (fun (name, raw_score) ->
      (name, sign *. (raw_score -. avg) /. std)
    ) name_raw_scores

let populate_ht spec =
  let norm_scores = normalize_scores spec in
  let ht = Ht.create (L.length norm_scores) in
  L.iter (fun (name, curr_score) ->
      try
        let prev_score = Ht.find ht name in
        if prev_score <> curr_score then
          begin
            Log.warn "%s: prev: %f curr: %f" name prev_score curr_score;
            let new_score =
              (* conflict resolution policy: keep only best score *)
              (if spec.increasing then min else max) prev_score curr_score in
            Ht.replace ht name new_score
          end
      with Not_found ->
        Ht.add ht name curr_score
    ) norm_scores;
  ht

let string_set_of_keys ht =
  StringSet.of_list (L.map fst (Ht.to_list ht))

type merge_policy = Min | Avg | Max

let merge_hts policy l =
  let reducer = match policy with
    | Avg -> BatList.favg
    | Min -> (BatList.min ~cmp:BatFloat.compare)
    | Max -> (BatList.max ~cmp:BatFloat.compare) in
  let all_names = L.map string_set_of_keys l in
  let common_names = match all_names with
    | [] -> failwith "Merge.merge_hts: no names"
    | x :: xs ->
      L.fold_left (fun acc y ->
          StringSet.inter acc y
        ) x xs in
  Log.info "common_names: %d" (StringSet.cardinal common_names);
  let name_avg_scores =
    StringSet.fold (fun name acc ->
        let scores = L.map (fun ht -> Ht.find ht name) l in
        (name, reducer scores) :: acc
      ) common_names [] in
  (* decreasing sort: we want high scores first *)
  L.sort (fun (_n1, s1) (_n2, s2) -> BatFloat.compare s2 s1) name_avg_scores

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -ifs <filename:name_field:{-}score_field>,...: \
              input score files\n  \
              the optional '-' means lower scores are better \
              (like ranks or docking scores)\n  \
              (field indexes start at 1)\n  \
              --avg: average score (default merging policy)\n  \
              --min: minimum score (merging policy)\n  \
              --max: maximum score (merging policy)\n  \
              -o <filename>: output scores file\n  \
              -n <int>: keep only top N\n  \
              -d <char>: field separator (default=\\t)\n"
       Sys.argv.(0);
     exit 1);
  let input_spec = CLI.get_string ["-ifs"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_top = CLI.get_int_opt ["-n"] args in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let policy = match (CLI.get_set_bool ["--min"] args,
                      CLI.get_set_bool ["--avg"] args,
                      CLI.get_set_bool ["--max"] args) with
  | (false, false, false) -> Avg (* default policy *)
  | (true, false, false) -> Min
  | (false, true, false) -> Avg
  | (false, false, true) -> Max
  | _ -> failwith "incompatible options: --min, --avg, --max" in
  CLI.finalize();
  let specs = parse_ifs sep input_spec in
  LO.with_out_file output_fn (fun out ->
      let hts = L.map populate_ht specs in
      let rescored = merge_hts policy hts in
      let to_write = match maybe_top with
        | None -> rescored
        | Some n -> Utls.list_really_take n rescored in
      L.iter (fun (name, score) ->
          Printf.fprintf out "%s\t%f\n" name score
        ) to_write
    )

let () = main ()
