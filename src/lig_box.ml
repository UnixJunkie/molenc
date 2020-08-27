(* given one ligand in a MOL2 file (xtal reference),
   compute geometric information about it
   (center, radius, bouding boxes, etc) *)

open Printf

module CLI = Minicli.CLI
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls
module V3 = Vector3

let index_of x l = match L.index_of x l with
  | Some i -> i
  | None -> assert(false)

let drop_pat_and_before pat l =
  let i = index_of pat l in
  L.drop (i + 1) l

let drop_pat_and_after pat l =
  let i = index_of pat l in
  L.take i l

let count pat l =
  L.length (L.filter ((=) pat) l)

let xyz_of_mol2_atom_line l =
  try
    (* " 1 O41 10.5170 3.9030 17.5470 O.2 1 O6K 0.0000" *)
    Scanf.sscanf l " %d %s %f %f %f %s %d %s %f"
      (fun _i _id x y z _type _j _res _charge -> x, y, z)
  with exn ->
    (Log.fatal "Lig_box.xyz_of_mol2_atom_line: cannot parse: %s" l;
     raise exn)

(* FBR: TODO output all info graphically in the BILD format for chimera *)

let fst3 (a, _, _) = a
let snd3 (_, b, _) = b
let trd3 (_, _, c) = c

let center xyz_l =
  let n = float (L.length xyz_l) in
  let init = V3.make 0.0 0.0 0.0 in
  let sum =
    L.fold_left (fun acc xyz ->
        V3.(add acc (of_triplet xyz))
      ) init xyz_l in
  V3.(to_triplet (div sum n))

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i ligand.mol2\n  \
              [-m <float>]: margin (in A) for all dimensions (default=0)\n"
       Sys.argv.(0);
     exit 1);
  Log.set_log_level Log.INFO;
  Log.set_output stderr;
  Log.color_on ();
  let input_fn = CLI.get_string ["-i"] args in
  let _margin = CLI.get_float_def ["-m"] args 0.0 in
  CLI.finalize ();
  let all_lines = Utls.lines_of_file input_fn in
  (* cout #molecules; enforce only 1 *)
  Utls.enforce (1 = count "@<TRIPOS>MOLECULE" all_lines)
    ("Lig_box: more than one molecule in " ^ input_fn);
  (* locate "@<TRIPOS>ATOM" line; delete all lines before *)
  let atom_lines' = drop_pat_and_before "@<TRIPOS>ATOM" all_lines in
  (* locate "@<TRIPOS>BOND" line; delete all lines after *)
  let atom_lines = drop_pat_and_after "@<TRIPOS>BOND" atom_lines' in
  let nb_atoms = L.length atom_lines in
  Log.info "%d atoms" nb_atoms;
  (* compute stuffs *)
  let triplets = L.map xyz_of_mol2_atom_line atom_lines in
  let xs, ys, zs = L.map fst3 triplets,
                   L.map snd3 triplets,
                   L.map trd3 triplets in
  let xmin, xmax = L.min_max xs in
  let ymin, ymax = L.min_max ys in
  let zmin, zmax = L.min_max zs in
  Log.info "x_{min|max}: %.2f %.2f" xmin xmax;
  Log.info "y_{min|max}: %.2f %.2f" ymin ymax;
  Log.info "z_{min|max}: %.2f %.2f" zmin zmax;
  let c_x, c_y, c_z = center triplets in
  Log.info "center: %.2f %.2f %.2f" c_x c_y c_z;
  let dx, dy, dz = xmax -. xmin, ymax -. ymin, zmax -. zmin in
  Log.info "dx dy dz: %.2f %.2f %.2f" dx dy dz;
  Log.info "englobing volume: %.2f A^3" (dx *. dy *. dz);
  ()

let () = main ()
