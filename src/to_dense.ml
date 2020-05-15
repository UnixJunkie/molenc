
(* read a molenc output file (.txt) and output it in dense csv format for R *)

open Printf

module CLI = Minicli.CLI
module Log = Dolog.Log
module Fp = Molenc.Fingerprint
module Utls = Molenc.Utls

let buff = Buffer.create 1024

let string_of_line nb_features line =
  Buffer.clear buff;
  try
    Scanf.sscanf line "%s@,%f,%s"
      (fun _name pIC50 fp_str ->
         bprintf buff "%f" pIC50;
         let fp = Fp.of_string fp_str in
         Fp.to_dense_bprintf buff nb_features fp;
         Buffer.contents buff
      )
  with exn ->
    (Log.error "cannot parse: %s" line;
     raise exn)

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  Log.info "start";
  let argc, args = CLI.init () in
  let show_help = CLI.get_set_bool ["-h";"--help"] args in
  if argc = 1 || show_help then
    (eprintf "usage:\n\
              %s [-np <int>] -i <molecules.txt> -n <nb_features>\n"
       Sys.argv.(0);
     exit 1);
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let input_fn = CLI.get_string ["-i"] args in
  let nb_features = CLI.get_int ["-n"] args in
  CLI.finalize ();
  (* header made of column numbers: IC50 in first column then features *)
  for i = 0 to nb_features do
    if i = 0 then printf "0"
    else printf " %d" i
  done;
  printf "\n";
  (* dense data lines *)
  Utls.iter_on_lines_of_file input_fn (fun line ->
      let str = string_of_line nb_features line in
      printf "%s\n" str
    )

let () = main ()
