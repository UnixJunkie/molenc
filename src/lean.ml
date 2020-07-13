
(* find the size of the random sample from the training set that is enough
   to train a model (can help avoiding overfitting, since we might not
   need to use the whole training set to tune the model.
   However, note that for production you will still need to train the
   model on the full training set.

   Bibliography:
   =============
   Domingos, P. (2012). "A few useful things to know about machine learning."
   Communications of the ACM, 55(10), 78-87.
 *)

open Printf

module CLI = Minicli.CLI
module Log = Dolog.Log
module F = BatFloat
module L = BatList
module Rand = BatRandom

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: input file\n  \
               -s <char>: field separator (default=\\t)\n  \
               -f <int>: field number (starts from 1)\n"
        Sys.argv.(0);
      exit 1
    end;
  let _sep = CLI.get_char_def ["-s"] args '\t' in
  let _field = CLI.get_int ["-f"] args in
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
  failwith "not implemented yet"

let () = main ()
