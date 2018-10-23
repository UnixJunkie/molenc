(* atom environment *)

open Printf

module L = MyList

(*   layer = (depth, counted-atoms) *)
type layer = int * ((PiEltHA.t * int) list)
type t = layer list

let counted_types_to_string (l: (PiEltHA.t * int) list): string =
  let buff = Buffer.create 80 in
  L.iteri (fun i (x, count) ->
      bprintf buff (if i = 0 then "%s:%d" else ",%s:%d")
        (PiEltHA.to_string x) count
    ) l;
  Buffer.contents buff

let counted_types_of_string (s: string): (PiEltHA.t * int) list =
  let strings = BatString.nsplit s ~by:"," in
  L.map (fun str -> Scanf.sscanf str "%s:%d" Utls.make_pair) strings

let layer_to_string ((depth, counted_types): layer): string =
  sprintf "%d-%s" depth (counted_types_to_string counted_types)

let layer_of_string (str: string): layer =
  Scanf.sscanf str "%d-%s" (fun d s ->
      (d, counted_types_of_string s)
    )

let to_string (layers: t): string =
  let buff = Buffer.create 80 in
  L.iteri (fun i layer ->
      bprintf buff (if i = 0 then "%s" else ";%s")
        (layer_to_string layer)
    ) layers;
  Buffer.contents buff

let of_string (s: string): t =
  let layer_strings = BatString.nsplit s ~by:";" in
  L.map layer_of_string layer_strings

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
