
module A = BatArray
module L = BatList
module Log = Dolog.Log

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_molecules l =
    create 1 Two_bands (A.of_list l)
end

(* For each molecule, find its nearest neighbor name and distance,
   over all Bsts *)
let nearest_neighbor_names ncores bst_fns mols =
  match bst_fns with
  | [] -> assert(false) (* at least one bst fn is required *)
  | fn :: fns ->
    let annot_mols =
      (* load one bst *)
      Log.info "loading %s..." fn;
      let (bst: Bstree.t) = Utls.restore fn in
      Parany.Parmap.parmap ncores
        (fun mol ->
           let nn, dist = Bstree.nearest_neighbor mol bst in
           (mol, FpMol.get_name nn, dist)
        ) mols in
    (* fold on the other BSTs *)
    L.fold_left (fun annotated bst_fn ->
        (* load another bst *)
        Log.info "loading %s..." bst_fn;
        let (bst: Bstree.t) = Utls.restore bst_fn in
        Parany.Parmap.parmap ncores (fun (mol, nn_name, dist) ->
            if dist = 0.0 then
              (* already nearest *)
              (mol, nn_name, dist)
            else
              let curr_nn, curr_dist = Bstree.nearest_neighbor mol bst in
              if curr_dist < dist then
                (mol, FpMol.get_name curr_nn, curr_dist)
              else
                (mol, nn_name, dist)
          ) annotated
      ) annot_mols fns
