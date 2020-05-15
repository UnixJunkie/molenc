(* Copyright (C) 2020, Francois Berenger

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
module Log = Dolog.Log
module Mini_mol = Molenc.Mini_mol
module Ht = BatHashtbl
module Scale = Molenc.Scale
module StringSet = BatSet.String
module Utls = Molenc.Utls

let read_one counter input () =
  try
    let m = Ap_types.read_one counter input in
    if !counter mod 1000 = 0 then
      (* user feedback *)
      eprintf "read %d\r%!" !counter;
    m
  with End_of_file ->
    begin
      Log.info "read %d" !counter;
      raise Parany.End_of_input
    end

let process_one radii m =
  let buff = Buffer.create 1024 in
  let name = Mini_mol.get_name m in
  bprintf buff "#%s\n" name;
  let seen_envs = Ht.create 1000 in
  L.iter (fun radius ->
      let envs = Mini_mol.encode radius m in
      L.iter (fun (env, count) ->
          (* only output envs that were not already encountered
             at lower radius *)
          if not (Ht.mem seen_envs env) then
            begin
              bprintf buff "%s %d\n" (Atom_env.to_string env) count;
              Ht.add seen_envs env ()
            end
        ) envs
    ) radii;
  Buffer.contents buff

let write_one output str =
  fprintf output "%s" str

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s -i molecules.{types|ph4} -r {radius|srad:frad} -o out.idx\n  \
              -i <filename>: where to read molecules from\n  \
              -r {<int>|<int>:<int>}: encoding radius or radii range\n  \
              -d <filename>: read feature dico from file\n  \
              -o <filename>: where to write encoded molecules\n  \
              [-n <int>]: max jobs in parallel\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-n"] args 1 in
  let scale =
    if L.mem "-r" args && L.mem "-d" args then
      (* enforce that radius ranges are equal *)
      let r_scale = Scale.of_string (CLI.get_string ["-r"] args) in
      let d_scale = Scale.of_dictionary_header (CLI.get_string ["-d"] args) in
      Utls.enforce (r_scale = d_scale)
        (sprintf "Encoder: -r and -d don't agree: r_scale=%s d_scale=%s"
           (Scale.to_string r_scale) (Scale.to_string d_scale));
      r_scale
    else
      match CLI.get_string_opt ["-r"] args with
      | Some r_str -> Scale.of_string r_str
      | None ->
        let dico_fn = CLI.get_string ["-d"] args in
        Scale.of_dictionary_header dico_fn in
  let radii = Scale.to_list scale in
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      (* format header *)
      fprintf output "#radius=%s\n%!" (Scale.to_string scale);
      Parany.run nprocs ~csize:1
        ~demux:(read_one (ref 0) input)
        ~work:(process_one radii)
        ~mux:(write_one output)
    )

let () = main ()
