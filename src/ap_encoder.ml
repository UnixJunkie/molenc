(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Atom pairs encoder *)

open Printf

type mode = Read_dictionary of string
          | Create_dictionary of string

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i molecules.types -o molecules.pairs\n  \
              -i <filename>: where to read molecules from\n  \
              -o <filename>: where to write encoded molecules\n  \
              -od <filename>: create and write encoding dictionary\n  \
                              (incompatible with -id)\n  \
              -id <filename>: read existing feature dictionary from file\n  
                              (incompatible with -od)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  failwith "not implemented yet"

let () = main ()
