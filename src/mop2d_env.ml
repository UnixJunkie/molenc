(* MOLPRINT 2D atom environment *)

open Printf

module L = MyList

(* center, neighbors 1 bond away, neighbors 2 bonds away, etc.
   MUST BE CANONICAL (i.e. sorted) *)
type t = Sybyl.t * (Sybyl.t * int) list list

let to_string ((center, neighbors): t): string =
  let counted_neighbors_str l =
    L.to_string (fun (typ, count) ->
        sprintf "(%s,%d)" (Sybyl.to_string typ) count
      ) l in
  sprintf "%s-%s"
    (Sybyl.to_string center)
    (L.to_string counted_neighbors_str neighbors)

let of_string (s: string): t =
  let center, neighbors =
    (* an atom env index has the count of how many times that atom env.
       was seen (upon creation of the index) but we ignore it *)
    try Scanf.sscanf s "%s@-%s@ %d" (fun a b _count -> a, b)
    with _ -> failwith ("Mop2d_env.of_string: cannot parse triplet: " ^ s) in
  let of_pair_str s' =
    try Scanf.sscanf s' "(%s@,%d)" (fun a b -> (Sybyl.of_string a, b))
    with _ -> failwith
                (sprintf "Mop2d_env.of_string: cannot parse pair: %s in %s"
                   s' s) in
  let neighbors' = BatString.chop ~l:1 ~r:1 neighbors in
  let sub_lists = BatString.nsplit neighbors' ~by:"];" in
  let sub_lists =
    L.map (fun s ->
        if BatString.ends_with s "]" then s
        else s ^ "]"
      ) sub_lists in
  (Sybyl.of_string center,
   L.map (L.of_string of_pair_str) sub_lists)

(* parse the 1st line of a .idx file *)
let parse_index_comment fn =
  let header, index_lines = Utls.maybe_extract_comment_header fn in
  match header with
  | None -> (-1, [])
  | Some comment ->
    let radius = Scanf.sscanf comment "#radius=%d" (fun r -> r) in
    (radius, index_lines)

(* parse the 1st line of a .mop2d file *)
let parse_molecules_comment fn =
  let header, mol_lines = Utls.maybe_extract_comment_header fn in
  match header with
  | None -> (-1, "/dev/null", mol_lines)
  | Some comment ->
    let radius, index_fn =
      Scanf.sscanf comment "#radius=%d;index=%s"
        (fun r fn -> (r, fn)) in
    (radius, index_fn, mol_lines)

(* parse the 1st line of an already opened .mop2d file
   (and advance the file pointer) *)
let parse_comment input =
  try (* we are parsing a valid .mop2d file *)
    Scanf.sscanf (input_line input)
      "#radius=%d;index=%s" (fun r fn -> (r, fn))
  with (* we are not *)
    Scanf.Scan_failure _ -> (-1, "/dev/null")

(* extract the MOP2D atom env. to bitstring index HT *)
let restore_mop2d_index fn =
  let radius, index_lines = parse_index_comment fn in
  let mop2d_envs = L.map of_string index_lines in
  let res = Hashtbl.create 11 in
  L.iteri (fun i env ->
      (* eprintf "%s\n" (Mop2d_env.to_string env); *)
      assert(not (Hashtbl.mem res env));
      Hashtbl.add res env i
    ) mop2d_envs;
  Log.info "index size: %d" (Hashtbl.length res);
  (radius, res)
