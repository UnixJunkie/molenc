(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* molecular encoder: a molecule is a list of atom environments.
   canonicalization is done by sorting (atoms in an environment
   as well as the list of environments that constitutes
   the entirely encoded molecule) *)

open Printf

module Ap_types = Molenc.Ap_types
module Atom_env = Molenc.Atom_env
module CLI = Minicli.CLI
module L = BatList
module Mini_mol = Molenc.Mini_mol
module Ht = BatHashtbl
module Scale = Molenc.Scale
module StringSet = BatSet.String
module Utls = Molenc.Utls

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.{types|ph4} -r {radius|srad:frad} -o output.idx\n\
              -i <filename>: where to read molecules from\n\
              -r {<int>|<int>:<int>}: encoding radius or radii range\n\
              -o <filename>: where to write encoded molecules\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  (* assert(BatString.ends_with input_fn ".types" ||
   *        BatString.ends_with input_fn ".ph4"); *)
  let output_fn = CLI.get_string ["-o"] args in
  let scale = Scale.of_string (CLI.get_string ["-r"] args) in
  let radii = Scale.to_list scale in
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      fprintf output "#radius=%s\n" (Scale.to_string scale);
      let counter = ref 0 in
      try
        while true do
          let m = Ap_types.read_one counter input in
          if !counter mod 1000 = 0 then
            eprintf "molecules seen: %d\r%!" !counter; (* user feedback *)
          let name = Mini_mol.get_name m in
          let seen_envs = Ht.create 11 in
          fprintf output "#%s\n" name;
          L.iter (fun radius ->
              let envs = Mini_mol.encode radius m in
              L.iter (fun (env, count) ->
                  (* only output envs that were not already encountered
                     at lower radius *)
                  if not (Ht.mem seen_envs env) then
                    (fprintf output "%s %d\n" (Atom_env.to_string env) count;
                     Ht.add seen_envs env ())
                ) envs
            ) radii
        done
      with End_of_file ->
        Log.info "read %d from %s" !counter input_fn
    )

let () = main ()
