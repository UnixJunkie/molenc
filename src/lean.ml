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
module Fn = Filename
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

(* [lst] must be sorted in increasing order
   [mini] must be <= min(lst)
   [maxi] must be >= max(lst) *)
let histo mini maxi n_steps lst =
  let bins = L.frange mini `To maxi n_steps in
  let avg x y =
    0.5 *. (x +. y) in
  let rec loop acc l = function
    | [] -> assert(false)
    | [_] -> acc
    | x :: y :: zs ->
      (* L.span is a kind of fold_while *)
      let this_bin, rest = L.span (fun x -> x < y) l in
      let n = L.length this_bin in
      loop ((avg x y, n) :: acc) rest (y :: zs) in
  let raw_histo = loop [] lst bins in
  (* normalize it; the Y-axis should be a frequency *)
  let total = float (L.length lst) in
  L.rev_map (fun (x, y) -> (x, (float y) /. total)) raw_histo

let histograms n_steps a1 a2 =
  let min1, max1 = A.min_max a1 in
  let min2, max2 = A.min_max a2 in
  let mini, maxi = (min min1 min2, max max1 max2) in
  let histo1 = histo mini maxi n_steps (A.to_list a1) in
  let histo2 = histo mini maxi n_steps (A.to_list a2) in
  (histo1, histo2)

let plot_histograms n_steps a1 a2 =
  let tmp_data_fn = Fn.temp_file "lean_histo_" ".txt" in
  Utls.with_out_file tmp_data_fn (fun out ->
      (* compute histograms *)
      let h1, h2 = histograms n_steps a1 a2 in
      (* dump them to a tmp file *)
      L.iter2 (fun (x1, y1) (x2, y2) ->
          assert(x1 = x2);
          fprintf out "%g %g %g\n" x1 y1 y2
        ) h1 h2
    );
  Log.info "tmp histo data: %s" tmp_data_fn;
  (* FBR:TODO call the venerable and venerated gnuplot *)
  ()

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
    let smaller_sample_size = A.length smaller_sample in
    total := smaller_sample_size + batch_size;
    let bigger_sample =
      A.of_list
        (L.map (get_field_as_float sep field)
           (L.take !total
              (L.shuffle ~state:(Rand.State.make_self_init ()) all_lines)))
    in
    A.sort BatFloat.compare smaller_sample;
    A.sort BatFloat.compare bigger_sample;
    plot_histograms 50 smaller_sample bigger_sample;
    let ks = Stats.ks2_test ~alpha smaller_sample bigger_sample in
    log_KS_test smaller_sample_size alpha ks
  done

let () = main ()
