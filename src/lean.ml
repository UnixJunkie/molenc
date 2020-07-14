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

module A = BatArray
module CLI = Minicli.CLI
module F = BatFloat
module L = BatList
module Log = Dolog.Log
module Rand = BatRandom
module Stats = Owl.Stats
module Utls = Molenc.Utls

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

let log_KS_test n alpha ks_res =
  let null_rejected = Owl.Stats.(ks_res.reject) in
  let p_value = Owl.Stats.(ks_res.p_value) in
  let ks_stat = Owl.Stats.(ks_res.score) in
  if p_value < alpha then
    (* distributions differ *)
    Log.warn "differ N=%d KS: %.3f p: %.3f (a=%.3f) rej=%d"
      n ks_stat p_value alpha (Utls.int_of_bool null_rejected)
  else
    (* distributions are the same *)
    Log.info "same N=%d KS: %.3f p: %.3f (a=%.3f) rej=%d"
      n ks_stat p_value alpha (Utls.int_of_bool null_rejected)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -b0 <int>: first batch size\n  \
               -b <int>: batch size (bootstrap size increment; \
               default=5%% of max)\n  \
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
  let all_lines = Utls.lines_of_file input_fn in
  let n = L.length all_lines in
  let batch_size =
    CLI.get_int_def ["-b"] args
      (int_of_float (F.ceil (0.05 *. (float n)))) in
  let batch_size0 = CLI.get_int_def ["-b0"] args batch_size in
  CLI.finalize ();
  Log.info "batch size: %d" batch_size;
  let total = ref 0 in
  while !total <= n do
    let smaller_sample =
      A.of_list
        (L.map (get_field_as_float sep field)
           (L.take (if !total = 0 then batch_size0 else !total)
              (L.shuffle ~state:(Rand.State.make_self_init ()) all_lines)))
    in
    total := (A.length smaller_sample) + batch_size;
    let bigger_sample =
      A.of_list
        (L.map (get_field_as_float sep field)
           (L.take !total
              (L.shuffle ~state:(Rand.State.make_self_init ()) all_lines)))
    in
    A.sort BatFloat.compare smaller_sample;
    A.sort BatFloat.compare bigger_sample;
    let ks = Stats.ks2_test ~alpha smaller_sample bigger_sample in
    log_KS_test !total alpha ks
  done

let () = main ()
