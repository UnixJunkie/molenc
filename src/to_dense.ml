
(* read a molenc output file (.txt) and output it in dense csv format for R *)

open Printf

module CLI = Minicli.CLI
module Log = Dolog.Log
module Fp = Molenc.Fingerprint
module Utls = Molenc.Utls

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  Log.info "start";
  let argc, args = CLI.init () in
  let show_help = CLI.get_set_bool ["-h";"--help"] args in
  if argc = 1 || show_help then
    (eprintf "usage:\n\
              %s -i <molecules.txt> -n <nb_features>\n"
       Sys.argv.(0);
     exit 1);
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
      try
        Scanf.sscanf line "%s@,%f,%s"
          (fun _name pIC50 fp_str ->
             printf "%f" pIC50;
             let a = Fp.to_dense nb_features (Fp.of_string fp_str) in
             for i = 0 to (Array.length a) - 1 do
               printf " %d" a.(i)
             done;
             printf "\n"
          )
      with exn ->
        (Log.error "cannot parse: %s" line;
         raise exn)
    )

let () = main ()
