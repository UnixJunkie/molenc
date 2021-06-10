(* Copyright (C) 2020, Francois Berenger

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
   Journal of Chemical Information and Computer Sciences, 39(4), 747-750. *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module BST = Bst.Bisec_tree.Make (FpMol)
module Ht = Hashtbl
module L = Molenc.MyList
module Log = Dolog.Log
module StringSet = BatSet.String
module Utls = Molenc.Utls

let sort_by_decr_density clusters =
  (* for repeatability of the clustering:
     stable_sort plus using cluster center names
     in case clusters have same size *)
  L.stable_sort (fun (c1, _m1, s1) (c2, _m2, s2) ->
      if s1 = s2 then
        compare (FpMol.get_name c1) (FpMol.get_name c2)
      else
        BatInt.compare s2 s1 (* decr sort *)
    ) clusters

let butina nprocs dist_t molecules =
  let rec loop clusters clustered = function
    | [] -> sort_by_decr_density clusters
    | ((center, members, size) :: rest) as remaining ->
      if size = 1 then
        (* all remaining clusters are singletons *)
        loop (L.rev_append remaining clusters) clustered []
      else
        let clusters' = (center, members, size) :: clusters in
        let clustered' = StringSet.union clustered members in
        (* update remaining clusters *)
        let to_cluster =
          L.filter (fun (c, _m, s) ->
              (* molecules previously included in a cluster
               * cannot be promoted to cluster center *)
              (s = 1 ||
               not (StringSet.mem (FpMol.get_name c) clustered'))
            ) rest in
        (* update their members and sizes *)
        let to_cluster' =
          L.map (fun (c, m, s) ->
              if s = 1 then (c, m, s)
              else
                let m' = StringSet.diff m clustered' in
                (c, m', StringSet.cardinal m')
            ) to_cluster in
        let to_cluster'' = sort_by_decr_density to_cluster' in
        loop clusters' clustered' to_cluster'' in
  (* density around each mol *)
  let bst = BST.(create 1 Two_bands (Array.of_list molecules)) in
  let mol_densities =
    Parany.Parmap.parmap ~csize:1 nprocs (fun mol ->
        let neighbors = BST.neighbors mol dist_t bst in
        let nb_neighbs = L.length neighbors in
        let neighbor_names = StringSet.of_list (L.map FpMol.get_name neighbors) in
        (mol, neighbor_names, nb_neighbs)
      ) molecules in
  let highest_density_first = sort_by_decr_density mol_densities in
  loop [] StringSet.empty highest_density_first

let list_for_all_vs_all p l =
  let rec loop = function
    | [] -> true
    | x :: xs ->
      L.for_all (p x) xs &&
      loop xs in
  loop l

(* check some mathematical properties about cluster centers and members *)
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
  (* no two cluster centers are nearer than d_t *)
  assert(list_for_all_vs_all (fun (c1, _m1, _s1) (c2, _m2, _s2) ->
      FpMol.dist c1 c2 >= d_t
    ) clusters)

(* FBR: do some sampling in order to find a reasonable Tanimoto threshold automatically *)

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
               [-t <float>]: Tanimoto threshold (default=0.9)\n  \
               [--verify]: run some checks on the obtained clusters\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let dist_t = 1.0 -. (CLI.get_float_def ["-t"] args 0.9) in
  assert(dist_t >= 0.0 && dist_t < 1.0);
  let verify = CLI.get_set_bool ["--verify"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize ();
  let all_mols = FpMol.molecules_of_file input_fn in
  let total = L.length all_mols in
  Log.info "read %d from %s" total input_fn;
  Log.info "BST: done";
  let clusters = butina nprocs dist_t all_mols in
  let single_count = L.filter_count (fun (_c, _m, s) -> s = 1) clusters in
  Log.info "clusters: %d" (L.length clusters);
  Log.info "singletons: %d" single_count;
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
      (* more checks *)
      if verify then
        verify_clusters dist_t all_mols clusters
    )

let () = main ()
