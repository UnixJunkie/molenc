(* Encoded molecules filtering and diversity selection.
   Functionalities:
   - remove training set molecules (or anything too near) from a
     database of molecules
     _OR_
   - enforce diversity in a database of molecules *)

open Printf

module CLI = Minicli.CLI

type mode = Filter of string (* file from where to read molecules to exclude *)
          | Diversify

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n \
               -i <filename>: molecules to filter (database)\n  \
               -o <filename>: output file\n  \
               [-t <float>]: Tanimoto threshold (default=1.0)\n  \
               [-e <filename>]: molecules to exclude from the ones to filter\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let threshold = CLI.get_float_def ["-t"] args 1.0 in
  let mode = match CLI.get_string_opt ["-e"] args with
    | None -> Diversify
    | Some fn -> Filter fn in
  CLI.finalize ();
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      let counter = ref 0 in
      try
        begin match mode with
        | Diversify -> assert(false)
        | Filter train_fn ->
          while true do
            failwith "not implemented yet"
          done
        end
      with End_of_file ->
        Log.info "read %d from %s" !counter input_fn
    )

let () = main ()
