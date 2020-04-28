(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Atom pairs encoder *)

open Printf

module Ap_types = Molenc.Ap_types
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module Log = Dolog.Log
module Mini_mol = Molenc.Mini_mol
module Utls = Molenc.Utls

type filename = string

type mode = Input_dictionary of filename
          | Output_dictionary of filename

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
        Scanf.sscanf line "%d %s" (fun id feat ->
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
    (if !counter mod 1000 = 0 then
       eprintf "read %d\r%!" !counter);
    m
  with End_of_file ->
    (Log.info "read %d" !counter;
     raise Parany.End_of_input)

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
    L.rev_map (fun (feat, count) ->
        (Ht.find feat2id feat, count)
      ) pairs in
  (* canonicalization *)
  let cano = L.sort compare_int_pairs feat_id_counts in
  L.iteri (fun i (feat_id, count) ->
      bprintf buff (if i > 0 then ";%d:%d" else "%d:%d")
        feat_id count
    ) cano;
  Buffer.add_char buff ']';
  Buffer.contents buff

let write_one output str =
  fprintf output "%s" str

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
              -np <int>: maximum number of cores\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let dico_mode = match (CLI.get_string_opt ["-id"] args,
                         CLI.get_string_opt ["-od"] args) with
  | (None, None) -> failwith "Ap_encoder: provide either -id or -od"
  | (Some _, Some _) -> failwith "Ap_encoder: -id and -od are exclusive"
  | (Some id_fn, None) -> Input_dictionary id_fn
  | (None, Some od_fn) -> Output_dictionary od_fn in
  CLI.finalize ();
  let dico = match dico_mode with
    | Input_dictionary id_fn -> dico_from_file id_fn
    | Output_dictionary od_fn -> failwith "not implemented yet" in
  failwith "not implemented yet"

let () = main ()
