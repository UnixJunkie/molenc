
(* Encode molecules from a dataset to k-bits.
   k is user-chosen.
   k distinct molecules from the dataset are used as vantage points and chosen
   randomly.

   TODO:
   - add --test option
     - draw 1000 molecule pairs randomly
     - output true Vs predicted Jaccard
   - assess estimated distance precision as a function of the number of bits k
   - compare performance to WMH
   - output a dictionary file at the end:
     - all vantage molecules
     - all thresholds

   Possible improvements:
   - make the encoding faster: there is not need for a bst now
     -> then check it goes faster
   - use Butina algo. to select vantage molecules?
   - use maximal diversity to select vantage molecules?
   - choose vantage molecules so that correlation with the true Jaccard
     is maximized?
   - choose thresholds so that correlation with the true Jaccard is maximized?
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

let in_count = ref 0

let read_one lref () = match !lref with
  | [] -> raise Parany.End_of_input
  | x :: xs ->
    let mol = FpMol.parse_one !in_count x in
    incr in_count;
    lref := xs;
    mol

let total = ref 0
let out_count = ref 0
let completed = ref 0

let write_one (name, bits) =
  printf "%s,0.0,%s\n" name (KBE.to_string bits);
  incr out_count;
  let current = (100 * !out_count) / !total in
  if current > !completed then
    begin
      eprintf "done: %d%%\r%!" current;
      completed := current
    end

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
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  (* read all molecules *)
  Log.info "reading molecules...";
  let all_lines = Utls.lines_of_file input_fn in
  total := L.length all_lines;
  Log.info "molecules: %d" !total;
  (* choose vantage molecules for encoding *)
  let vmols =
    let rand_k_lines =
      L.really_take k
        (L.shuffle ~state:(BatRandom.State.make_self_init ()) all_lines) in
    L.mapi FpMol.parse_one rand_k_lines in
  let ht, thresholds, bst = KBE.init vmols in
  let dt, () =
    Utls.time_it (fun () ->
        Parany.run ~verbose:false ~csize:100 ~nprocs
          ~demux:(read_one (ref all_lines))
          ~work:(KBE.encode ht thresholds bst)
          ~mux:write_one
      ) in
  let rate = (float !total) /. dt in
  Log.info "rate: %.2f mols/s" rate

let () = main ()
