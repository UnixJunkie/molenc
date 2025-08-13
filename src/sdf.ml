
module Utls = Molenc.Utls

(* one molecule in SDF format (i.e. consecutive lines from a .sdf file) *)
type t = string

let line_buffer = Buffer.create 1024

(* Gzip.input_line does not exist yet:
   https://github.com/xavierleroy/camlzip/issues/51
   WARNING: very inefficient
*)
let gzip_input_line (input: Gzip.in_channel): string =
  try
    let c = ref (Gzip.input_char input) in
    while !c <> '\n' do
      Buffer.add_char line_buffer !c;
      c := Gzip.input_char input
    done;
    let res = Buffer.contents line_buffer in
    Buffer.clear line_buffer;
    res
  with End_of_file ->
    let res = Buffer.contents line_buffer in
    Buffer.clear line_buffer;
    if res = "" then
      raise End_of_file
    else
      res

exception Read_one

let read_one (input: in_channel): t =
  let buff = Buffer.create 10240 in
  try
    while true do
      let line = input_line input in
      if line = "$$$$" then (* end of molecule in SDF format *)
        (Buffer.add_string buff line;
         Buffer.add_char buff '\n';
         raise Read_one)
      else
        (Buffer.add_string buff line;
         Buffer.add_char buff '\n')
    done;
    assert(false)
  with
  | End_of_file | Read_one ->
    let res = Buffer.contents buff in
    if res = "" then
      raise End_of_file
    else
      res

let read_one_zip (input: Gzip.in_channel): t =
  let buff = Buffer.create 10240 in
  try
    while true do
      let line = gzip_input_line input in
      if line = "$$$$" then (* end of molecule in SDF format *)
        (Buffer.add_string buff line;
         Buffer.add_char buff '\n';
         raise Read_one)
      else
        (Buffer.add_string buff line;
         Buffer.add_char buff '\n')
    done;
    assert(false)
  with
  | End_of_file | Read_one ->
    let res = Buffer.contents buff in
    if res = "" then
      raise End_of_file
    else
      res

(* return the inchi string, no trailing '\n' *)
let get_inchi (mol: t): string =
  let line_before = "> <PUBCHEM_IUPAC_INCHI>\n" in
  let n = String.length line_before in
  try
    let i = BatString.find mol line_before in
    let j = i + n in
    let k = BatString.find_from mol j "\n" in
    BatString.sub mol j (k - j)
  with Not_found ->
    failwith ("Sdf.get_inchi: no inchi for: " ^ mol)

let get_inchikey (mol: t): string =
  let line_before = "> <PUBCHEM_IUPAC_INCHIKEY>\n" in
  let n = String.length line_before in
  try
    let i = BatString.find mol line_before in
    let j = i + n in
    let k = BatString.find_from mol j "\n" in
    BatString.sub mol j (k - j)
  with Not_found ->
    failwith ("Sdf.get_inchikey: no inchikey for: " ^ mol)

let get_fst_line m =
  fst (BatString.split m ~by:"\n")
