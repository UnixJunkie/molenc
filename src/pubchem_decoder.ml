(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* decode pubchem FPs (881 bits) to liblinear format
   Here is a pubchem FP example line:
---
4036230,0.0,0001010110...01010100001
---
 *)

open Printf

module CLI = Minicli.CLI
module Log = Dolog.Log
module String = BatString
module Utls = Molenc.Utls

let liblinear_line_of_pubchem_line line =
  match String.split_on_char ',' line with
  | [name; _IC50; bitstring] ->
    let is_active = String.starts_with name "active" in
    let nb_bits = String.length bitstring in
    assert(nb_bits = 881 || nb_bits = 16384);
    let buff = Buffer.create 16500 in
    Buffer.add_string buff (if is_active then "+1" else "-1");
    String.iteri (fun i c ->
        if c = '1' then
          let k = i + 1 in (* in liblinear: feature indexes start at 1 *)
          Printf.bprintf buff " %d:1" k
      ) bitstring;
    Buffer.contents buff
  | _ -> failwith ("Pubchem_decoder: invalide line: " ^ line)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage: %s\n  \
              -i <filename>: encoded molecules\n  \
              -o <filename>: decoded molecules for liblinear\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  CLI.finalize ();
  let line_counter = ref 0 in
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      try
        while true do
          let in_line = input_line input in
          incr line_counter;
          if !line_counter mod 1000 = 0 then
            eprintf "read: %d\r%!" !line_counter;
          let out_line = liblinear_line_of_pubchem_line in_line in
          fprintf output "%s\n" out_line
        done
      with End_of_file -> ()
    );
  eprintf "read: %d\n" !line_counter

let () = main ()
