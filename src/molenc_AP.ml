(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, Tokyo University,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
 *
 * New (2023) counted atom pairs encoder w/ optional molecular standardization. *)

open Printf

module CLI = Minicli.CLI
module LO = Line_oriented
module Log = Dolog.Log

type mode = Input_dict of string
          | Output_dict of string

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  (if argc = 1 then
     (eprintf "usage:\n  \
               %s -i in.smi -o out.AP\n  \
               -i <input.smi>: input molecules\n  \
               -o <output.AP>: unfolded counted atom pairs output\n  \
               -d <dico.dix>: use existing feature dictionary\n  \
               [-m <int>]: maximum atom pairs distance (in bonds; default=OFF)\n  \
               [-np <int>]: parallelize on NCORES (default=1)\n  \
               [-cs <int>]: chunk size (default=50)\n"
        Sys.argv.(0);
      exit 1)
  );
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let _csize = CLI.get_int_def ["-cs"] args 50 in
  let _max_dist_opt = CLI.get_int_opt ["-m"] args in
  let _dico_mode =
    match CLI.get_string_opt ["-d"] args with
    | None -> Output_dict (input_fn ^ ".dix")
    | Some fn -> Input_dict fn in
  CLI.finalize (); (* --------------------------------------------------- *)
  LO.with_infile_outfile input_fn output_fn (fun _input _output ->
      Parany.run ~preserve:true nprocs
        ~demux:(fun _ -> failwith "not implemented yet")
        ~work:(fun _ -> failwith "not implemented yet")
        ~mux:(fun _ -> failwith "not implemented yet")
    )

let () = main ()
