
(* prune features from a dictionary file, provided a list of feature indexes
 * to remove *)

module CLI = Minicli.CLI
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module S = BatString
module Utls = Molenc.Utls

open Printf

let prune_dico features_to_drop in_dico_fn out_dico_fn =
  let n = L.length features_to_drop in
  Utls.enforce (n > 0) "Model.prune_dico: |features_to_drop| = 0";
  let to_drop = Ht.create n in
  L.iter (fun i ->
      Ht.add to_drop i ()
    ) features_to_drop;
  Log.info "pruning %d features" (Ht.length to_drop);
  Utls.with_out_file out_dico_fn (fun out ->
      Utls.iter_on_lines_of_file in_dico_fn (fun line ->
          if S.starts_with line "#" then
            (* preserve comments *)
            fprintf out "%s\n" line
          else
            try
              (* this is the format for an Atom Pairs dictionary *)
              Scanf.sscanf line "%d %s@ %d" (fun featId featStr featCount ->
                  (* drop some features *)
                  if not (Ht.mem to_drop featId) then
                    fprintf out "%d %s %d\n" featId featStr featCount
                )
            with exn -> (Log.fatal "dico %s: cannot parse line %s"
                           in_dico_fn line;
                         raise exn)
        )
    )

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i input.dix -o output.dix -f features.txt\n  \
              -i <filename>: input AP dictionary\n  \
              -o <filename>: output AP dictionary\n  \
              -f <filename>: file with list of features to drop\n  \
              (format: one integer per line)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-i"] args in
  let features_fn = CLI.get_string ["-f"] args in
  CLI.finalize();
  let features = Utls.map_on_lines_of_file features_fn int_of_string in
  prune_dico features input_fn output_fn
