(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Work with DeepSMILES *)

open Printf

module A = Array
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module RNG = BatRandom.State
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString

(* because the Rdkit module uses Pyml *)
let () = Py.initialize ~version:3 ()

type mol = Rdkit.t

let mol_of_smiles (smi: string): mol =
  Rdkit.__init__ ~smi ()

let tok0 = "PAD"
let tok1 = "START"
let tok2 = "UNK"
let unknown_index = 2

let large_int = 1073741823 (* (2^30) - 1 *)

let lookup_in_dico dico word =
  try Ht.find dico word
  with Not_found ->
    let new_idx = Ht.length dico in
    Ht.add dico word new_idx;
    new_idx

let update_word_count word2count word =
  let prev = Ht.find_default word2count word 0 in
  Ht.replace word2count word (prev + 1)

let insert_reserved_tokens dico =
  let idx0 = lookup_in_dico dico tok0 in
  assert(idx0 = 0);
  let idx1 = lookup_in_dico dico tok1 in
  assert(idx1 = 1);
  let idx2 = lookup_in_dico dico tok2 in
  assert(idx2 = 2)

(* most frequent words have lower indexes *)
let dico_of_word2count word2count =
  let len = Ht.length word2count in
  let n = len + 3 in (* there are reserved tokens *)
  let ht = Ht.create n in
  insert_reserved_tokens ht;
  let word_counts = Ht.bindings word2count in
  let sorted =
    (* decr. sort *)
    L.sort (fun (_w1, c1) (_w2, c2) -> BatInt.compare c2 c1) word_counts in
  L.iter (fun (word, _count) ->
      let _i: int = lookup_in_dico ht word in
      ()
    ) sorted;
  ht

let dico_to_file fn dico =
  Log.info "saving dico to %s" fn;
  LO.with_out_file fn (fun out ->
      let kvs = Ht.bindings dico in
      let sorted = L.sort (fun (_w1, i) (_w2, j) -> BatInt.compare i j) kvs in
      L.iter (fun (word, index) ->
          fprintf out "%d %s\n" index word
        ) sorted
    )

type encoding = Bitstring (* unfolded-uncounted fingerprint / bag of words *)
              | Count_vector (* unfolded-counted fingerprint *)
              | Sequence (* vector of word indexes *)

let encoding_of_string = function
  | "bits" -> Bitstring
  | "count" -> Count_vector
  | "seq" -> Sequence
  | other -> failwith ("Dsmi.encoding_of_string: not supported: " ^ other)

let bits_of_dsmi dico dsmi =
  let n = Ht.length dico in
  let toks = S.split_on_char ' ' dsmi in
  let bits = Bytes.init n (fun _ -> '0') in
  L.iter (fun tok ->
      let i = Ht.find_default dico tok unknown_index in
      Bytes.set bits i '1'
    ) toks;
  let s = Bytes.unsafe_to_string bits in
  let buff = Buffer.create 1024 in
  (* AP file format *)
  Buffer.add_char buff '[';
  let start = ref true in
  S.iteri (fun i c ->
      if c = '1' then
        (if !start then
           (start := false;
            Printf.bprintf buff "%d:1" i)
         else
           Printf.bprintf buff ";%d:1" i)
    ) s;
  Buffer.add_char buff ']';
  Buffer.contents buff

let counts_of_dsmi dico dsmi =
  let n = Ht.length dico in
  let toks = S.split_on_char ' ' dsmi in
  let counts = A.make n 0 in
  L.iter (fun tok ->
      let i = Ht.find_default dico tok unknown_index in
      counts.(i) <- counts.(i) + 1
    ) toks;
  let buff = Buffer.create 1024 in
  (* AP file format *)
  Buffer.add_char buff '[';
  let start = ref true in
  A.iteri (fun i c ->
      if c > 0 then
        (if !start then
           (start := false;
            Printf.bprintf buff "%d:%d" i c)
         else
           Printf.bprintf buff ";%d:%d" i c)
    ) counts;
  Buffer.add_char buff ']';
  Buffer.contents buff

let main () =
  let start = Unix.gettimeofday () in
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename>]: SMILES input file\n  \
              [-s <int>]: RNG seed\n  \
              [-n <int>]: DeepSMILES per input SMILES (default=1)\n  \
              -o <filename>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-e <string>]: output encoding (bits|count|seq)\n  \
              [--rand]: SMILES randomization (even if n=1)\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let maybe_seed = CLI.get_int_opt ["-s"] args in
  let n = CLI.get_int_def ["-n"] args 1 in
  let randomize = (CLI.get_set_bool ["--rand"] args) || n > 1 in
  let output_fn = CLI.get_string ["-o"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let encoding = CLI.get_string_def ["-e"] args "bits" in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let encoding_style = encoding_of_string encoding in
  let rng = match maybe_seed with
    | Some s -> RNG.make [|s|] (* repeatable *)
    | None -> RNG.make_self_init () in
  let nb_molecules = ref 0 in
  let dummy = mol_of_smiles "C" in
  let word_count = Ht.create 997 in (* enough to hold all DeepSMILES tokens *)
  let all_dsmi = ref [] in
  Log.info "reading SMILES...";
  LO.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let line = input_line input in
          incr nb_molecules;
          let smi, name = S.split line ~by:"\t" in
          let seed = RNG.int rng large_int in
          let dsmiles = Rdkit.get_deep_smiles dummy ~seed ~n ~randomize ~smi () in
          A.iter (Printf.fprintf output "%s\n") dsmiles;
          (* update word count *)
          A.iter (fun deep_smile ->
              all_dsmi := (name, deep_smile) :: !all_dsmi;
              let tokens = S.split_on_char ' ' deep_smile in
              L.iter (update_word_count word_count) tokens
            ) dsmiles
        done
      with End_of_file -> ()
    );
  all_dsmi := L.rev !all_dsmi; (* put them in input order *)
  Log.info "creating dictionary...";
  let dico = dico_of_word2count word_count in
  Log.info "seen %d words (3 are reserved)" (Ht.length dico);  
  let dico_fn = input_fn ^ ".dix" in
  dico_to_file dico_fn dico;
  Log.info "encoding...";
  let code_out_fn = output_fn ^ ".code" in
  (match encoding_style with
   | Bitstring ->
      LO.with_out_file code_out_fn (fun out ->
          L.iter (fun (name, dsmi) ->
              let bitstring = bits_of_dsmi dico dsmi in
              fprintf out "%s,0.0,%s\n" name bitstring
            ) !all_dsmi
        )
   | Count_vector ->
      LO.with_out_file code_out_fn (fun out ->
          L.iter (fun (name, dsmi) ->
              let feat_counts = counts_of_dsmi dico dsmi in
              fprintf out "%s,0.0,%s\n" name feat_counts
            ) !all_dsmi
        )
   | _ -> failwith "not implemented yet"
  );
  let stop = Unix.gettimeofday () in
  let dt = stop -. start in
  Log.info "encoding: %.1f molecule/s" ((float !nb_molecules) /. dt)

let () = main ()
