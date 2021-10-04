
module A = BatArray
module L = BatList
module Log = Dolog.Log

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_molecules l =
    create 1 Two_bands (A.of_list l)
end

(* For each molecule, find its nearest neighbor name and distance,
   over all Bsts; parallelized over molecules. *)
let nearest_neighbor_names ncores bst_fns mols =
  match bst_fns with
  | [] -> []
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

let bst_nearest_name_dist bst mol =
  let nn, dist = Bstree.nearest_neighbor mol bst in
  (FpMol.get_name nn, dist)

(* For each molecule, find its nearest neighbor name and distance,
   over all Bsts; parallelized over bisector trees (indexed chunks) *)
let nearest_neighbor_names_a
    (ncores: int) (bst_fns: string list) (mols_a: FpMol.t array)
  : (FpMol.t * string * float) array =
  match bst_fns with
  | [] -> [||]
  | fn :: fns' ->
    (* init accumulater, for the muxer process.
       This is the only calculation parallelized over molecules.
       Remaining calculations will be paralellized over Bsts. *)
    let annot_mols =
      Log.info "loading %s..." fn;
      let (bst: Bstree.t) = Utls.restore fn in
      Parany.Parmap.array_parmap ncores
        (fun mol ->
           let name, dist = bst_nearest_name_dist bst mol in
           (mol, name, dist))
        (mols_a.(0), "", 1.0)
        mols_a in
    let fns = A.of_list fns' in
    let () =
      Parany.run ncores
        ~demux:(
          let i = ref 1 in (* fn already processed *)
          let n = A.length fns in
          fun () ->
            if !i < n then
              let res = !i in
              incr i;
              res
            else
              raise Parany.End_of_input
        )
        ~work:(fun i ->
            Log.info "loading %s..." fns.(i);
            let (bst: Bstree.t) = Utls.restore fns.(i) in
            A.map (bst_nearest_name_dist bst) mols_a
          )
        ~mux:(
          let m = A.length mols_a in
          fun nearest_name_dists ->
            assert(A.length nearest_name_dists = m);
            for i = 0 to m - 1 do
              let mol, _prev_nearest_name, prev_dist =
                A.unsafe_get annot_mols i in
              if prev_dist = 0.0 then
                (* already nearest *) ()
              else
                let curr_nearest_name, curr_dist =
                  A.unsafe_get nearest_name_dists i in
                if curr_dist < prev_dist then (* update acc *)
                  A.unsafe_set annot_mols i (mol, curr_nearest_name, curr_dist)
            done
        ) in
    annot_mols
