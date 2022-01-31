(* Copyright (C) 2021, Francois Berenger

   Tsuda laboratory, The University of Tokyo, Japan.

   Bounding-Box Applicability Domain for Counted Atom Pairs
   (uses 2nd generation AP typing scheme) *)

open Printf

module A = Array
module CLI = Minicli.CLI
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module Rdkit = Molenc.Rdkit.Rdkit
module S = BatString

let () = Py.initialize () (* because the Rdkit module uses Pyml *)
(* FBR: try to embed the Python source file statically into the exe
   asked in pyml on github *)

(* (atomic_num, nb_HA, nb_H, HA_used_val, formal_charge) cf. rdkit_wrapper.py *)
type atom = (int * int * int * int * int)
type pair = (atom * int * atom)
type mol = Rdkit.t

(* sorting the types for canonicalization *)
let create_pair (src: atom) (dist: int) (dst: atom): pair =
  if compare src dst <= 0 then
    (src, dist, dst)
  else
    (dst, dist, src)

let string_of_atom ((a,b,c,d,e): atom): string =
  sprintf "(%d,%d,%d,%d,%d)" a b c d e

let atom_of_string (s: string): atom =
  try Scanf.sscanf s "(%d,%d,%d,%d,%d)" (fun a b c d e -> (a,b,c,d,e))
  with exn ->
    let () = Log.fatal "AP_BBAD.atom_of_string: cannot parse: %s" s in
    raise exn

let string_of_pair ((src, dist, dst): pair): string =
  Printf.sprintf "%s-%d-%s" (string_of_atom src) dist (string_of_atom dst)

let pair_of_string (s: string): pair =
  try Scanf.sscanf s "(%d,%d,%d,%d,%d)-%d-(%d,%d,%d,%d,%d)"
        (fun a b c d e dist f g h i j ->
           ((a,b,c,d,e), dist, (f,g,h,i,j)))
  with exn ->
    let () = Log.fatal "AP_BBAD.pair_of_string: cannot parse: %s" s in
    raise exn

let feature_count_from_line line =
  (* AP-BBAD line format:
     ^(6,2,2,2,0)-23-(6,3,0,4,0) 174$ *)
  try Scanf.sscanf line "%s@ %d"
        (fun feat max_count -> (pair_of_string feat, max_count))
  with exn ->
    let () = Log.fatal
        "AP_BBAD.feature_count_from_line: cannot parse: %s" line in
    raise exn

let bbad_from_file (fn: string): (pair, int) Ht.t =
  let lines = LO.lines_of_file fn in
  let n = L.length lines in
  let bbad = Ht.create n in
  L.iter (fun line ->
      let feat, max_count = feature_count_from_line line in
      Ht.add bbad feat max_count
    ) lines;
  bbad

let in_AP_BBAD bbad mol =
  Ht.for_all (fun feat count ->
      let max_count = Ht.find_default bbad feat 0 in
      count <= max_count
    ) mol

let mol_of_smiles (smi: string): mol =
  Rdkit.__init__ ~smi ()

let type_atom (mol: mol) (i: int): atom =
  let typ = Rdkit.type_atom mol ~i:i () in
  assert(A.length typ == 5);
  (typ.(0), typ.(1), typ.(2), typ.(3), typ.(4))

(* counted atom pairs *)
let encode_molecule (mol: mol): (pair, int) Ht.t =
  let n = Rdkit.get_num_atoms mol () in
  let ht = Ht.create (n * (n - 1) / 2) in
  for i = 0 to n - 1 do
    let src = type_atom mol i in
    for j = i to n - 1 do
      let dst = type_atom mol j in
      let dist = Rdkit.get_distance mol ~i ~j () in
      let feature = create_pair src dist dst in
      let prev_count = Ht.find_default ht feature 0 in
      Ht.replace ht feature (prev_count + 1)
    done
  done;
  ht

let read_one_line input () =
  try input_line input
  with End_of_file -> raise Parany.End_of_input

let process_one line =
  let smi, _name = S.split line ~by:"\t" in
  let mol = mol_of_smiles smi in
  encode_molecule mol

let annotate_one line =
  let smi, _name = S.split line ~by:"\t" in
  let mol = mol_of_smiles smi in
  (encode_molecule mol, line)

let filter_one bbad line =
  let mol = process_one line in
  if in_AP_BBAD bbad mol then
    Some line
  else
    None

let write_one total in_AD out = function
  | None -> incr total
  | Some line ->
    begin
      incr total;
      incr in_AD;
      fprintf out "%s\n" line
    end

let record_one count res pairs =
  incr count;
  res := pairs :: !res

(* Create a new AD using the union of two existing ones *)
let union_AD (ad1: (pair, int) Ht.t) (ad2: (pair, int) Ht.t)
  : (pair, int) Ht.t =
  Ht.merge (fun _feat maybe_count1 maybe_count2 ->
      match (maybe_count1, maybe_count2) with
      | (Some n, None  )
      | (None  , Some n) -> Some n
      | (Some m, Some n) -> Some (BatInt.max m n)
      | (None  , None  ) -> assert(false)
    ) ad1 ad2

let write_AD_out output_fn bbad =
  (* canonical form: write BBAD out decr. sorted by feature_count *)
  let key_values = A.of_list (Ht.to_list bbad) in
  A.sort (fun (key1, n1) (key2, n2) ->
      let cmp = BatInt.compare n2 n1 in
      if cmp = 0 then
        compare key1 key2 (* define a total order *)
      else
          cmp
    ) key_values;
  LO.with_out_file output_fn (fun out ->
      A.iter (fun (feat, max_count) ->
          fprintf out "%s %d\n" (string_of_pair feat) max_count
        ) key_values
    )

let main () =
  Log.color_on ();
  Log.(set_log_level INFO);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              [-i <filename>]: SMILES input file\n  \
              -o <filename>: output file\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [--bbad <fn1[,fn2[,...]]>]: previously computed BBAD\n  \
              to apply as filter\n  \
              will compute the union of ADs if several filenames\n  \
              [-v]: verbose/debug mode\n" Sys.argv.(0);
     exit 1);
  let verbose = CLI.get_set_bool ["-v"] args in
  if verbose then Log.(set_log_level DEBUG);
  let input_fn = CLI.get_string_def ["-i"] args "/dev/stdin" in
  let output_fn = CLI.get_string ["-o"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let bbad_fn = CLI.get_string_opt ["--bbad"] args in
  CLI.finalize ();
  match bbad_fn with
  | Some fn ->
    if S.contains fn ',' then
      let fns = S.split_on_string ~by:"," fn in
      match fns with
      | [] -> assert(false)
      | fn1 :: others ->
        let ad1 = bbad_from_file fn1 in
        let union_AD =
          L.fold (fun acc fn' ->
              let ad = bbad_from_file fn' in
              union_AD acc ad
            ) ad1 others in
        write_AD_out output_fn union_AD
    else
      (* if BBAD fn provided: apply BBAD on set of molecules
         - any unknown feature -> molecule filtered out
         - feature value too high -> molecule filtered out *)
      let bbad = bbad_from_file fn in
      let nb_molecules = ref 0 in
      let in_AD = ref 0 in
      let start = Unix.gettimeofday () in
      LO.with_in_file input_fn (fun input ->
          LO.with_out_file output_fn (fun out ->
              Parany.run ~preserve:true nprocs
                ~demux:(read_one_line input)
                ~work:(filter_one bbad)
                ~mux:(write_one nb_molecules in_AD out);
            )
        );
      let stop = Unix.gettimeofday () in
      let dt = stop -. start in
      Log.info "in_AD: %d/%d" !in_AD !nb_molecules;
      Log.info "filtering: %.1f molecule/s" ((float !nb_molecules) /. dt)
  | None ->
    (* default mode: compute BBAD for a set of molecules *)
    let nb_molecules = ref 0 in
    let start = Unix.gettimeofday () in
    let atom_pairs = ref [] in
    LO.with_in_file input_fn (fun input ->
        Parany.run nprocs ~demux:(read_one_line input)
          ~work:process_one
          ~mux:(record_one nb_molecules atom_pairs)
      );
    let stop = Unix.gettimeofday () in
    let dt = stop -. start in
    Log.info "encoding: %.1f molecule/s" ((float !nb_molecules) /. dt);
    let bbad = Ht.create !nb_molecules in
    L.iter (
      Ht.iter (fun feat count ->
          let prev_count = Ht.find_default bbad feat count in
          Ht.replace bbad feat (max count prev_count)
        )
    ) !atom_pairs;
    write_AD_out output_fn bbad

let () = main ()
