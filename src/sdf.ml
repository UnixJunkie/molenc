
(* one molecule in SDF format (i.e. consecutive lines from a .sdf file) *)
type t = string

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
