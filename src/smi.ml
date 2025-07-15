
(* read one molecule from a SMILES file *)
let read_one (input: in_channel): string =
  let line = input_line input in
  (* strip protects against trailing '\r'
     we append '\n' because all other formats
     end molecules with a '\n' *)
  (BatString.strip line) ^ "\n"

let read_one_zip (_: Gzip.in_channel) =
  failwith "Smi.read_one_zip: not implemented yet"

let get_name smiles_line =
  let _smiles, name =
    try BatString.split smiles_line ~by:"\t"
    with Not_found -> failwith "Smi.get_name: smiles file not using tabs" in
  (* all other file formats expect molecule name only,
   * not molecule name followed by EOL *)
  if BatString.ends_with name "\r\n" then
    BatString.rchop ~n:2 name
  else if BatString.ends_with name "\n" then
    BatString.rchop ~n:1 name
  else
    name
