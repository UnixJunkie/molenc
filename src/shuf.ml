(* replacement for UNIX's shuf command, but seedable for reproducibility *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module LO = Line_oriented
module Log = Dolog.Log
module RNG = BatRandom.State

exception Early_stop

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file\n  \
               -o <filename>: output file\n  \
               [-n <int>]: output at most N lines (default=all)\n  \
               [-s <int>]: random seed (default=none)\n"
        Sys.argv.(0);
      exit 1
    end;
  CLI.finalize (); (* ------------------------------------------------------ *)
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_n = CLI.get_int_opt ["-n"] args in
  let maybe_seed = CLI.get_int_opt ["-s"] args in
  (* input *)
  let all_lines = A.of_list (LO.lines_of_file input_fn) in
  let count = A.length all_lines in
  (* output all or not? *)
  let output_n = match maybe_n with
    | None -> count
    | Some m -> min count m in
  let rng = match maybe_seed with
    | None -> RNG.make_self_init ()
    | Some s -> RNG.make [|s|] in
  (* shuffle *)
  A.shuffle ~state:rng all_lines;
  (* output *)
  LO.with_out_file output_fn (fun out ->
      try
        A.iteri (fun i line ->
            if i < output_n then
              fprintf out "%s\n" line
            else
              raise Early_stop
          ) all_lines
      with Early_stop -> ()
    )

let () = main ()
