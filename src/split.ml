(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   Split large files with many molecules into chunks of given size.
   .smi, .sdf and .mol2 are supported. *)

open Printf

module A = Array
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

type file_format = SMILES
                 | SDF
                 | Mol2

let file_format_of_filename fn =
  if S.ends_with fn ".smi" then
    (SMILES, ".smi")
  else if S.ends_with fn ".smiles" then
    (SMILES, ".smiles")
  else if S.ends_with fn ".sdf" then
    (SDF, ".sdf")
  else if S.ends_with fn ".mol2" then
    (Mol2, ".mol2")
  else
    failwith
      ("Split.file_format_of_filename: not {.smi|.smiles|.sdf|.mol2}: %s" ^ fn)

(* molecule separators *)
let mol2_header = "@<TRIPOS>MOLECULE"
let sdf_footer = "$$$$"

exception Read_one

let read_one ff count input =
  match ff with
  | Mol2 ->
    begin
      let first_line = input_line input in
      incr count;
      let res =
        if first_line = mol2_header then
          ref [mol2_header]
        else
          (* add back mol2 header from previous molecule read *)
          ref [first_line; mol2_header] in
      try
        while true do
          let line = input_line input in
          if line = mol2_header then
            raise Read_one (* consumed the mol2 header from next molecule *)
          else
            res := line :: !res
        done;
        assert(false)
      with End_of_file | Read_one -> L.rev !res
    end
  | SMILES ->
    begin
      try
        let line = input_line input in
        incr count;
        [line]
      with End_of_file -> []
    end
  | SDF ->
    let res = ref [] in
    try
      while true do
        let line = input_line input in
        if line = sdf_footer then
          begin
            res := line :: !res;
            incr count;
            raise Read_one
          end
        else
          res := line :: !res;
      done;
      assert(false)
    with Read_one | End_of_file -> L.rev !res

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -i <filename.{smi|mol2|sdf}>: molecules input file\n  \
              [-c <int>]: chunk size (molecules per output file; default=50)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let csize = CLI.get_int_def ["-c"] args 50 in
  CLI.finalize ();
  let ff, ext = file_format_of_filename input_fn in
  let basename = Filename.chop_suffix input_fn ext in
  let count = ref 0 in
  LO.with_in_file input_fn (fun input ->
      let output_fn = sprintf "%s_%04d%s" basename !count ext in
      LO.with_out_file output_fn (fun out ->
          try
            for _i = 1 to csize do
              let mol = read_one ff count input in
              L.iter (fprintf out "%s\n") mol
            done
          with End_of_file -> ()
        )
    );
  Log.info "read %d from %s" !count input_fn

let () = main ()
