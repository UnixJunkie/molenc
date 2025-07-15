(* Copyright (C) 2022, Francois Berenger

   Tsuda Laboratory, Graduate School of Frontier Sciences,
   The University of Tokyo, Japan.

   Support .ph4 files. *)

(* A .ph4 files has format:
---
<num_feats:int>:<mol_name:string>
ARO 1.47088 -0.706617 1.86095
.
.
.
--- *)

module S = BatString

let parse_header_line line =
  let num_feats, name = S.split ~by:":" line in
  (int_of_string num_feats, name)

exception Read_one

(* read one molecule from a .ph4 file *)
let read_one (input: in_channel): string =
  let buff = Buffer.create 2048 in
  try
    let line = input_line input in
    let num_feats, _name = parse_header_line line in
    Buffer.add_string buff line;
    Buffer.add_char buff '\n';
    for _i = 1 to num_feats do
      let line = input_line input in
      Buffer.add_string buff line;
      Buffer.add_char buff '\n'
    done;
    raise Read_one
  with | End_of_file | Read_one ->
    let res = Buffer.contents buff in
    if res = "" then
      raise End_of_file
    else res


let read_one_zip (_: Gzip.in_channel) =
  failwith "Ph4.read_one_zip: not implemented yet"

let get_name ph4_lines =
  let header, _rest = S.split ph4_lines ~by:"\n" in
  let _num_feats, name = parse_header_line header in
  name
