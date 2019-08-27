(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* A fingerprint-encoded molecule *)

module A = Array
module Fp = Fingerprint
module Ht = Hashtbl
module L = MyList

type t = { name: string;
           index: int; (* position in input file *)
           fp: Fp.t }

let create name index bitstring =
  { name; index; fp = Fp.of_string bitstring }

(* read one molecule from an FP file *)
let read_one_mol line =
  try Scanf.sscanf line "%s@,%f,%s"
        (fun name _ic50 bitstring ->
           (name, bitstring)
        )
  with Scanf.Scan_failure msg ->
    failwith ("FpMol.read_one_mol: fmt: %s@,%f,%s err: " ^ msg ^
              " line: " ^ line)

let parse_one index line =
  let name, bitstring = read_one_mol line in
  create name index bitstring

let molecules_of_file fn =
  Utls.mapi_on_lines_of_file fn parse_one

let dist m1 m2 =
  Fp.distance m1.fp m2.fp

let get_name x =
  x.name

let get_index x =
  x.index

let get_fp x =
  x.fp

let max_feature_index x =
  Fp.max_feature_index x.fp

let mol_is_active line =
  BatString.starts_with line "active"

let is_active x =
  mol_is_active x.name
