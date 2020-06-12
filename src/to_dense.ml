
(* read a molenc output file (.txt) and output it in dense csv format for R *)

open Printf

module A = BatArray
module Bloom = Molenc.Bloom
module CLI = Minicli.CLI
module Log = Dolog.Log
module Fp = Molenc.Fingerprint
module Utls = Molenc.Utls

let expand_line nb_features maybe_bloom line =
  try
    Scanf.sscanf line "%s@,%f,%s"
      (fun _name pIC50 fp_str ->
         printf "%f" pIC50;
         let fp = Fp.of_string fp_str in
         (match maybe_bloom with
          | None -> Fp.to_dense_printf nb_features fp;
          | Some params ->
            let bloom = Bloom.encode params fp in
            A.iter (printf " %d") bloom
         );
         printf "\n"
      )
  with exn ->
    (Log.fatal "cannot parse: %s" line;
     raise exn)

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  Log.info "start";
  let argc, args = CLI.init () in
  let show_help = CLI.get_set_bool ["-h";"--help"] args in
  if argc = 1 || show_help then
    (eprintf "usage:\n\
              %s [-np <int>] -i <molecules.txt> -n <nb_features>\n\
              [--bloom <int>,<int>]: k,m counted Bloom filter params\n"
       Sys.argv.(0);
     exit 1);
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let input_fn = CLI.get_string ["-i"] args in
  let input_features = CLI.get_int ["-n"] args in
  let output_features, maybe_bloom =
    match CLI.get_string_opt ["--bloom"] args with
    | None -> (input_features, None)
    | Some k_m ->
      Scanf.sscanf k_m "%d,%d" (fun k m ->
          Utls.enforce (m < input_features) "m >= input_features";
          (m, Some (Bloom.init input_features k m))
        ) in
  CLI.finalize ();
  (* CSV header made of column numbers: IC50 in first column then features *)
  printf "0";
  for i = 1 to output_features do
    printf " %d" i
  done;
  printf "\n";
  (* dense data lines, with optional counted Bloom filter encoding *)
  Utls.iter_on_lines_of_file input_fn
    (expand_line input_features maybe_bloom)

let () = main ()
