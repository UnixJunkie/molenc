
(* initial tests of random coordinates sub-sampling wihout replacement *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = BatList
module Utls = Molenc.Utls

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i encoded_molecules.txt\n" Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  (* read all molecules *)
  let molecules = FpMol.molecules_of_file input_fn in
  let nb_features = L.max (L.map FpMol.nb_features molecules) in
  let fingerprints = A.of_list (L.map FpMol.get_fp molecules) in
  let drop_rates = L.frange 0.1 `To 0.9 9 in
  let feat_id_max = nb_features - 1 in
  let all_feat_ids = L.range 0 `To feat_id_max in
  L.iter (fun drop_p ->
      failwith "not implemented yet"
    ) drop_rates

let () = main ()
