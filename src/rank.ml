(* compute ranks associated to each score
   equal scores will be given equal ranks *)

open Printf

module CLI = Minicli.CLI
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module LO = Line_oriented
module String = BatString
module Utls = Molenc.Utls

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input scores file\n  \
               -o <filename>: output rank and scores file\n  \
               -f <int>: score field (>= 1)\n  \
               -d <char>: field separator (default=\\t)\n \
               [-r]: increasing scores order (default=decreasing)\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let sep = CLI.get_char_def ["-d"] args '\t' in
  let field = (CLI.get_int ["-f"] args) - 1 in
  let revert = CLI.get_set_bool ["-r"] args in
  CLI.finalize();
  let all_scores = ref [] in
  (* read all scores *)
  Utls.iter_on_lines_of_file input_fn (fun line ->
      let score_field = String.cut_on_char sep field line in
      let score =
        try Scanf.sscanf score_field "%f" (fun x -> x)
        with exn ->
          begin
            Log.fatal "Rank: cannot parse float: %s" score_field;
            raise exn
          end in
      all_scores := score :: !all_scores
    );
  (* create the score to rank LUT *)
  let uniq_scores =
    let cmp =
      if revert then BatFloat.compare (* increasing sort *)
      else
        (* default: scores in decreasing order; i.e. the highest score
           gets the lowest rank *)
        (fun x y -> BatFloat.compare y x) in
    L.sort_uniq cmp !all_scores in
  let score2rank = Ht.create (L.length uniq_scores) in
  L.iteri (fun i score ->
      Ht.add score2rank score i
    ) uniq_scores;
  (* output all lines, allong with their rank *)
  LO.with_out_file output_fn (fun output ->
      Utls.iter_on_lines_of_file input_fn (fun line ->
          let score_field = String.cut_on_char sep field line in
          let score =
            try Scanf.sscanf score_field "%f" (fun x -> x)
            with exn ->
              begin
                Log.fatal "Rank: cannot parse float: %s" score_field;
                raise exn
              end in
          let rank = Ht.find score2rank score in
          fprintf output "%s%c%d\n" line sep rank
        )
    )

let () = main ()
