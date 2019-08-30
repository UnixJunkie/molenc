
(* Encode molecules from a dataset to k-bits.
   k is user-chosen.
   k distinct molecules from the dataset are used as vantage points and chosen
   randomly.

   TODO:
   - parallelize encoder using Parany
   - assess estimated distance precision as a function of k
   - compare performance to WMH

   Possible improvements:
   - choose thresholds randomly?
   - use Butina algo. to select vantage molecules?
   - use maximal diversity to select vantage molecules?
   - choose vantage molecules so that correlation with the true Jaccard
     is maximized?
 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Fp = Molenc.Fingerprint
module FpMol = Molenc.FpMol
module L = Molenc.MyList
module Utls = Molenc.Utls
module WMH = Molenc.WMH
module KBE = Molenc.KBE

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i encoded_molecules.txt -k nbits\n" Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let k = CLI.get_int ["-k"] args in
  (* read all molecules *)
  Log.info "reading molecules...";
  let all_molecules = FpMol.molecules_of_file input_fn in
  let n = L.length all_molecules in
  Log.info "read %d" n;
  (* choose vantage molecules for encoding *)
  let vmols =
    L.really_take k
      (L.shuffle ~state:(BatRandom.State.make_self_init ()) all_molecules) in
  let ht, bst = KBE.init vmols in
  L.iter (fun mol ->
      let name = FpMol.get_name mol in
      let kbe = KBE.encode ht bst mol in
      printf "%s,0.0,%s\n" name (KBE.to_string kbe)
    ) all_molecules

let () = main ()
