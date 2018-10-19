
module At = Atom
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

(* throw End_of_file once there is no more to read *)
let read_one (counter: int ref) (input: in_channel): Mini_mol.t =
  let consume_until tag input =
    let rec loop () =
      let line = input_line input in
      if tag = line then ()
      else loop ()
    in
    loop ()
  in
  (* skip lines until we start to read a new molecule *)
  consume_until molecule_header input;
  (* the line just after is the molecule name *)
  let mol_name = input_line input in
  (* the line after tells how many atoms and bonds there are *)
  let nb_atoms, nb_bonds = read_header (input_line input) in
  let graph = Array.make nb_atoms Node.dummy in
  (* skip lines until we start to read atoms *)
  consume_until atoms_header input;
  (* read all atoms *)
  for _ = 1 to nb_atoms do
    let a = At.of_mol2_line (input_line input) in
    graph.(At.(a.idx)) <- Node.create At.(a.typ) IntSet.empty
  done;
  let line = input_line input in
  assert(line = bonds_header);
  (* read all bonds *)
  for _ = 1 to nb_bonds do
    let b = Bond.of_mol2_line (input_line input) in
    let src = Bond.(b.src) in
    let dst = Bond.(b.dst) in
    let curr = graph.(src) in
    let curr' = graph.(dst) in
    graph.(src) <- Node.add_succ curr dst;
    (* edges are directed in a MOL2 file;
       while the real molecular graph is undirected *)
    graph.(dst) <- Node.add_succ curr' src
  done;
  let res = Mini_mol.create mol_name graph in
  incr counter;
  res

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
