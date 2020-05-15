(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Atom pairs encoder *)

open Printf

module Ap_types = Molenc.Ap_types
module Atom_pair = Molenc.Atom_pair
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module Log = Dolog.Log
module Mini_mol = Molenc.Mini_mol
module Utls = Molenc.Utls

type mode = Input_dictionary of string
          | Output_dictionary of string

(* reconstruct the map feat->featId from given file *)
let dico_from_file fn =
  let n = Utls.count_lines_of_file fn in
  assert(n > 1);
  let feat2id = Ht.create (n - 1) in
  let header = Utls.get_first_line fn in
  Utls.enforce (header = "#atom_pairs")
    ("Ap_encoder.dico_from_file: not an atom pairs dict: " ^ fn);
  Utls.iter_on_lines_of_file fn (fun line ->
      if not (BatString.starts_with line "#") then
        Scanf.sscanf line "%d %s %d" (fun id feat _count ->
            (* the binding defined in the dictionary should be unique *)
            assert(not (Ht.mem feat2id feat));
            Ht.add feat2id feat id
          )
    );
  feat2id

let read_one counter input () =
  try
    let m = Ap_types.read_one counter input in
    (* user feedback *)
    (if !counter mod 1000 = 0
     then eprintf "read %d\r%!" !counter);
    m
  with End_of_file ->
    (Log.info "read %d" !counter;
     raise Parany.End_of_input)

(* create and store the feature dictionary *)
let dico_to_file molecules_fn dico_fn =
  (* How many times a given feature is seen when considering all
   * molecules (note that each molecule may have it several times) *)
  let feat_to_tot_count = Ht.create 100_000 in
  Utls.with_in_file molecules_fn (fun input ->
      try
        let mol_counter = ref 0 in
        while true do
          let mol = read_one mol_counter input () in
          let pair_counts = Mini_mol.atom_pairs mol in
          L.iter (fun (pair, curr_count) ->
              let feat_str = Atom_pair.to_string pair in
              let prev_count = Ht.find_default feat_to_tot_count feat_str 0 in
              Ht.replace feat_to_tot_count feat_str (prev_count + curr_count)
            ) pair_counts
        done
      with Parany.End_of_input -> ()
    );
  (* make the feature dictionary invariant to input molecules order
     by sorting (i.e. canonicalization) *)
  let sorted =
    let count_feats =
      let feat_counts = Ht.bindings feat_to_tot_count in
      L.rev_map (fun (feat_str, count) -> (count, feat_str)) feat_counts in
    Utls.list_rev_sort compare count_feats in
  let n = L.length sorted in
  let dico = Ht.create n in
  Utls.with_out_file dico_fn (fun output ->
      fprintf output "#atom_pairs\n";
      L.iteri (fun id (count, feat_str) ->
          (* we also print the total count, to allow verification
           * of the dictionary's features order *)
          fprintf output "%d %s %d\n" id feat_str count;
          assert(not (Ht.mem dico feat_str)); (* feature must be uniq *)
          Ht.add dico feat_str id
        ) sorted
    );
  dico

(* specialized compare for int pairs *)
let compare_int_pairs (i, j) (k, l) =
  let cmp = BatInt.compare i k in
  if cmp <> 0 then cmp
  else BatInt.compare j l

let process_one feat2id mol =
  let buff = Buffer.create 1024 in
  let name = Mini_mol.get_name mol in
  bprintf buff "%s,0.0,[" name;
  let pairs = Mini_mol.atom_pairs mol in
  let feat_id_counts =
    L.fold_left (fun acc (feat, count) ->
        let feat_str = Atom_pair.to_string feat in
        try (Ht.find feat2id feat_str, count) :: acc
        with Not_found ->
          (Log.warn "Ap_encoder: mol: %s unknown feat: %s %d"
             name feat_str count;
           acc)
      ) [] pairs in
  (* canonicalization *)
  let cano = L.sort compare_int_pairs feat_id_counts in
  L.iteri (fun i (feat_id, count) ->
      bprintf buff (if i > 0 then ";%d:%d" else "%d:%d")
        feat_id count
    ) cano;
  Buffer.add_char buff ']';
  Buffer.contents buff

let write_one output str =
  fprintf output "%s\n" str

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i molecules.types -o molecules.pairs {-od|-id} ap.dix\n  \
              -i <filename>: input molecules file\n  \
              -o <filename>: encoded molecules output file\n  \
              -od <filename>: create and write feature dictionary\n      \
              (incompatible with -id)\n  \
              -id <filename>: read existing feature dictionary\n      \
              (incompatible with -od)\n  \
              -np <int>: max number of cores (default=1)\n  \
              -cs <int>: chunk size (default=1)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-cs"] args 1 in
  let dico_mode =
    match (CLI.get_string_opt ["-id"] args,
           CLI.get_string_opt ["-od"] args) with
    | (None, None) -> failwith "Ap_encoder: provide either -id or -od"
    | (Some _, Some _) -> failwith "Ap_encoder: -id and -od are exclusive"
    | (Some id_fn, None) -> Input_dictionary id_fn
    | (None, Some od_fn) -> Output_dictionary od_fn in
  CLI.finalize ();
  let dico = match dico_mode with
    | Input_dictionary id_fn -> dico_from_file id_fn
    | Output_dictionary od_fn -> dico_to_file input_fn od_fn in
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      Parany.run nprocs ~csize
        ~demux:(read_one (ref 0) input)
        ~work:(process_one dico)
        ~mux:(write_one output)
    )

let () = main ()
