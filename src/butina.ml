(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   Clustering using the Butina algorithm.

   Butina, D. (1999).
   "Unsupervised data base clustering based on daylight's
   fingerprint and Tanimoto similarity: A fast and automated way to cluster
   small and large data sets".
   Journal of Chemical Information and Computer Sciences, 39(4), 747-750.

   -> output cluster centers
   -> output cluster membership for each molecule *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module BST = Bst.Bisec_tree.Make (FpMol)
module Ht = Hashtbl
module L = BatList
module StringSet = BatSet.String
module Utls = Molenc.Utls

let verbose = ref false

let sort_by_decr_density clusters =
  L.sort (fun (_c1, _m1, s1) (_c2, _m2, s2) ->
      BatInt.compare s2 s1
    ) clusters

let butina nprocs dist_t molecules =
  let rec loop clusters clustered = function
    | [] -> L.rev clusters
    | (center, members, size) :: rest ->
      let clusters' = (center, members, size) :: clusters in
      let clustered' = StringSet.union clustered members in
      (* update remaining clusters *)
      let to_cluster =
        L.filter (fun (c, _m, _s) ->
            not (StringSet.mem (FpMol.get_name c) clustered')
          ) rest in
      (* update their members and sizes *)
      let to_cluster' =
        L.map (fun (c, m, _s) ->
            let m' = StringSet.diff m clustered' in
            (c, m', StringSet.cardinal m')
          ) to_cluster in
      let to_cluster'' = sort_by_decr_density to_cluster' in
      loop clusters' clustered' to_cluster'' in
  (* density around each mol *)
  let bst = BST.(create 1 Two_bands (Array.of_list molecules)) in
  let mol_densities =
    Parmap.parmap ~ncores:nprocs ~chunksize:1 (fun mol ->
        let neighbors = BST.neighbors mol dist_t bst in
        let nb_neighbs = L.length neighbors in
        let neighbor_names = StringSet.of_list (L.map FpMol.get_name neighbors) in
        (mol, neighbor_names, nb_neighbs)
      ) (Parmap.L molecules) in
  (* isolated molecules (singleton clusters) are handled here *)
  let singletons, others = L.partition (fun (_c, _m, s) -> s = 1) mol_densities in
  let clustered =
    List.fold_left (fun acc (c, _m, _s) ->
        StringSet.add (FpMol.get_name c) acc
      ) StringSet.empty singletons in
  (* update members *)
  let to_cluster =
    L.map (fun (c, m, _s) ->
        let m' = StringSet.diff m clustered in
        (c, m', StringSet.cardinal m')
      ) others in
  let highest_density_first = sort_by_decr_density to_cluster in
  let clusters = loop [] StringSet.empty highest_density_first in
  L.append clusters singletons

let verify_clusters d_t molecules clusters =
  let n = L.length molecules in
  let name2mol = Ht.create n in
  L.iter (fun m ->
      Ht.add name2mol (FpMol.get_name m) m
    ) molecules;
  (* all molecules from one cluster are at d <= d_t from its center *)
  assert(L.for_all (fun (center_mol, member_names, _size) ->
      StringSet.for_all (fun mol_name ->
          let mol = Ht.find name2mol mol_name in
          FpMol.dist center_mol mol <= d_t
        ) member_names
    ) clusters);
  (* no two cluster centers are nearer than 2 * d_t *)
  (* TODO *)
  ()

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s\n  \
               -i <filename>: molecules to filter (\"database\")\n  \
               -o <filename>: output file\n  \
               [-t <float>]: Tanimoto threshold (default=1.0)\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let dist_t = 1.0 -. (CLI.get_float_def ["-t"] args 0.9) in
  assert(dist_t >= 0.0 && dist_t < 1.0);
  verbose := CLI.get_set_bool ["-v"] args;
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize ();
  let all_mols = FpMol.molecules_of_file input_fn in
  let total = L.length all_mols in
  Log.info "read %d from %s" total input_fn;
  Log.info "BST OK";
  let clusters = butina nprocs dist_t all_mols in
  Log.info "clusters: %d" (L.length clusters);
  Utls.with_out_file output_fn (fun out ->
      let clustered =
        L.fold_lefti (fun count cid (center, members, size) ->
            fprintf out "cid: %d size: %d center: %s members:"
              cid size (FpMol.get_name center);
            StringSet.iter (fprintf out " %s") members;
            fprintf out "\n";
            (count + size)
          ) 0 clusters in
      (* ensure all molecules were clustered *)
      assert(clustered = total);
      (* more checks FBR: disable those tests by default *)
      verify_clusters dist_t all_mols clusters
    )

let () = main ()
