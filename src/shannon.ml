(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Compute Shannon information entropy of atom pair features *)

open Printf

module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module Ht = BatHashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

let log2_scale = 1.0 /. (log 2.0)

let log2 x =
  log2_scale *. (log x)

let shannon_entropy n val_counts =
  let res = ref 0.0 in
  Ht.iter (fun _value count ->
      let p_i = (float count) /. n in
      res := !res +. (p_i *. log2 p_i)
    ) val_counts;
  -1.0 *. !res

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i mols.txt -n num_features -o entropy.txt\n  \
              -i <filename>: input molecules file\n  \
              -n <int>: number of features\n  \
              -o <filename>: output file\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let nb_features = CLI.get_int ["-n"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize ();
  Log.info "reading molecules...";
  let all_molecules =
    Utls.map_on_lines_of_file input_fn (FpMol.parse_one 0) in
  Log.info "read: %d" (L.length all_molecules);
  Log.info "computing entropy...";
  let ht = Ht.create nb_features in
  for i = 0 to nb_features - 1 do
    Ht.add ht i (Ht.create 11)
  done;
  L.iter (fun mol ->
      Fp.kv_iter (fun feat_id feat_count ->
          let acc = Ht.find ht feat_id in
          let prev_count = Ht.find_default acc feat_count 0 in
          Ht.replace acc feat_count (prev_count + 1)
        ) (FpMol.get_fp mol)
    ) all_molecules;
  let total = ref 0 in
  for i = 0 to nb_features - 1 do
    let acc = Ht.find ht i in
    Ht.iter (fun _k v ->
        total := !total + v
      ) acc
  done;
  Log.info "total: %d" !total;
  let n = float !total in
  (* entropy of each feature *)
  let feat_ent = ref [] in
  for i = 0 to nb_features - 1 do
    let acc = Ht.find ht i in
    let ent = shannon_entropy n acc in
    if ent > 0.0 then
      feat_ent := (i, ent) :: !feat_ent
  done;
  (* sort features by decreasing entropy *)
  let feat_encr_decr =
    L.sort (fun (_i, ei) (_j, ej) ->
        BatFloat.compare ej ei
      ) !feat_ent in
  let cumulated = ref 0.0 in
  Utls.with_out_file output_fn (fun out ->
      L.iter (fun (feat, ent) ->
          fprintf out "%d %f %f\n" feat ent !cumulated;
          cumulated := !cumulated +. ent
        ) feat_encr_decr
    )

let () = main ()
