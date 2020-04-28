(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Atom pairs encoder *)

open Printf

module CLI = Minicli.CLI
module Ht = BatHashtbl
module Log = Dolog.Log
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
              (incompatible with -od)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
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
