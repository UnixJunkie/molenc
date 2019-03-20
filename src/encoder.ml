
(* output atom environments; sorted by decreasing frequency given
   a file containing molecules *)

open Printf

module CLI = Minicli.CLI
module L = MyList
module Ht = BatHashtbl
module StringSet = BatSet.String

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.{types|ph4} -r max_radius -o output.idx\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  assert(BatString.ends_with input_fn ".types" ||
         BatString.ends_with input_fn ".ph4");
  let output_fn = CLI.get_string ["-o"] args in
  let radius = CLI.get_int ["-r"] args in
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      fprintf output "#radius=%d\n" radius;
      let counter = ref 0 in
      try
        while true do
          let m = Ap_types.read_one counter input in
          if !counter mod 1000 = 0 then
            eprintf "%d molecules seen\r%!" !counter; (* user feedback *)
          let envs = Mini_mol.encode radius m in
          L.iter (fun (env, count) ->
              fprintf output "%s %d\n" (Atom_env.to_string env) count
            ) envs
        done
      with End_of_file ->
        Log.info "read %d from %s" !counter input_fn
    )

let () = main ()
