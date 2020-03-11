(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

module At = Syb_atom
module IntSet = BatSet.Int
module S = BatString

let molecule_header = "@<TRIPOS>MOLECULE"
let atoms_header = "@<TRIPOS>ATOM"
let bonds_header = "@<TRIPOS>BOND"

(* parse line just after the molecule name line *)
let read_header l =
  try Scanf.sscanf l " %d %d %d %d %d "
        (fun nb_atoms nb_bonds _ _ _ -> (nb_atoms, nb_bonds))
  with _ -> failwith ("read_header: could not parse: " ^ l)

exception Read_one

let buff = Buffer.create 10240

(* read one molecule from a MOL2 file *)
let read_one_raw (input: in_channel): string =
  try
    while true do
      let line = input_line input in
      if line = molecule_header && Buffer.length buff <> 0 then
        (* just finished reading one *)
        raise Read_one
      else
        (Buffer.add_string buff line;
         Buffer.add_char buff '\n')
    done;
    assert(false)
  with
  | Read_one ->
    let res = Buffer.contents buff in
    Buffer.reset buff;
    Buffer.add_string buff molecule_header; (* put in buffer next mol's header *)
    Buffer.add_char buff '\n';
    res
  | End_of_file ->
    if Buffer.length buff = 0 then
      raise End_of_file
    else
      let res = Buffer.contents buff in
      Buffer.reset buff;
      res

let get_name mol_lines =
  let _header, rest = S.split mol_lines ~by:"\n" in
  let name, _tail = S.split rest ~by:"\n" in
  name
