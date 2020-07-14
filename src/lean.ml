(*
   find the size of the random sample from the training set that is enough
   to train a model (can help avoiding overfitting, since we might not
   need to use the whole training set to tune the model).
   However, note that for production you will still need to train your
   model on the full training set.

   Bibliography:
   =============
   Domingos, P. (2012). "A few useful things to know about machine learning."
   Communications of the ACM, 55(10), 78-87.
 *)

open Printf

module A = Array
module CLI = Minicli.CLI
module F = BatFloat
module L = BatList
module Log = Dolog.Log
module Rand = BatRandom
module Stats = Owl.Stats

let get_field_as_float s f line =
  try
    let field_str = L.at (BatString.split_on_char s line) f in
    Scanf.sscanf field_str "%f" (fun x -> x)
  with exn ->
    let () =
      Log.fatal "Lean.get_field_as_float: \
                 cannot parse float from field %d with sep '%c' \
                 in line: %s" f s line in
    raise exn

let log_KS_test n Owl.Stats.{reject; p_value; score} =
  if reject then
    Log.warn "N=%d dists differ.; KS: %.3f p-value: %.3f"
      n p_value score
  else
    Log.info "N=%d same dist.; KS: %.3f p-value: %.3f"
      n p_value score

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -b <int>: batch size (default=5%% of max)\n  \
               -i <filename>: input file\n  \
               -s <char>: field separator (default=\\t)\n  \
               -f <int>: field number (starts from 1)\n  \
               -a <float>: KS-test alpha (default=0.05)\n"
        Sys.argv.(0);
      exit 1
    end;
  let alpha = CLI.get_float_def ["-a"] args 0.05 in
  let sep = CLI.get_char_def ["-s"] args '\t' in
  let field = (CLI.get_int ["-f"] args) - 1 in
  let input_fn = CLI.get_string ["-i"] args in
  let rand_lines =
    let all_lines = Molenc.Utls.lines_of_file input_fn in
    L.shuffle ~state:(Rand.State.make_self_init ()) all_lines in
  let n = L.length rand_lines in
  let batch_size =
    CLI.get_int_def ["-b"] args
      (int_of_float (F.ceil (0.05 *. (float n)))) in
  CLI.finalize ();
  Log.info "batch size: %d" batch_size;
  let all_values =
    A.of_list
      (L.map (get_field_as_float sep field) rand_lines) in
  A.sort BatFloat.compare all_values; (* sort in increasing order *)
  let curr_sample = ref (L.take batch_size rand_lines) in
  let next_sample = ref (L.take (2 * batch_size) rand_lines) in
  let total = ref (L.length !curr_sample) in
  while !total <= n do
    let smaller_sample = 
      A.of_list (L.map (get_field_as_float sep field) !curr_sample) in
    let bigger_sample =
      A.of_list (L.map (get_field_as_float sep field) !next_sample) in
    A.sort BatFloat.compare smaller_sample;
    A.sort BatFloat.compare bigger_sample;
    let ks = Stats.ks2_test ~alpha:alpha smaller_sample all_values in
    log_KS_test !total ks;
    (* let ks = Stats.ks2_test ~alpha:alpha smaller_sample bigger_sample in
     * log_KS_test !total ks; *)
    total := !total + batch_size;
    curr_sample := !next_sample;
    next_sample := L.take !total rand_lines
  done

let () = main ()
