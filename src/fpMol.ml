(* Copyright (C) 2020, Francois Berenger

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

open Printf

type t = { name: string;
           index: int; (* position in input file *)
           value: float;
           fp: Fp.t }

let create name index value bitstring =
  { name; index; value; fp = Fp.of_string bitstring }

(* read one molecule from an FP file *)
let read_one_mol line =
  try Scanf.sscanf line "%s@,%f,%s"
        (fun name value bitstring ->
           (name, value, bitstring)
        )
  with Scanf.Scan_failure msg ->
    failwith ("FpMol.read_one_mol: fmt: %s@,%f,%s err: " ^ msg ^
              " line: " ^ line)

let parse_one index line =
  let name, value, bitstring = read_one_mol line in
  create name index value bitstring

(* go back to the line format you came from *)
let to_string (m: t): string =
  sprintf "%s,%g,[%s]"
    m.name
    m.value
    (Fp.to_string m.fp)

let to_out out m =
  fprintf out "%s\n" (to_string m)

let molecules_of_file fn =
  Utls.mapi_on_lines_of_file fn parse_one

let dist m1 m2 =
  Fp.distance m1.fp m2.fp

let tani m1 m2 =
  Fp.tanimoto m1.fp m2.fp

let get_name x =
  x.name

let get_value x =
  x.value

let get_index x =
  x.index

let get_fp x =
  x.fp

let nb_features x =
  Fp.nb_features x.fp

let mol_is_active line =
  BatString.starts_with line "active"

let is_active x =
  mol_is_active x.name

let drop_features to_drop x =
  { x with fp = Fp.drop_features to_drop x.fp }
