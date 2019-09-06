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
module Utls = Molenc.Utls

let verbose = ref false

let butina tani_t molecules =
  let n = L.length molecules in
  let clusters = Ht.create n in
  let clustered = Ht.create n in
  let dist_t = 1.0 -. tani_t in
  let cid = ref 0 in
  let rec loop mols =
    match mols with
    | [] -> ()
    | _ ->
      (* assign neighborhood density to each mol *)
      let bst = BST.(create 20 One_band (Array.of_list mols)) in
      (* FBR: detect and remove all singletons here *)
      let mol_densities =
        (* FBR: this should use parmap *)
        L.map (fun mol ->
            let neighbors = BST.neighbors mol dist_t bst in
            (mol, neighbors, L.length neighbors)
          ) mols in
      let highest_density_first =
        L.sort (fun (_m1, _mn1, n1) (_m2, _mn2, n2) ->
            BatInt.compare n2 n1
          ) mol_densities in
      match highest_density_first with
      | [] -> assert(false)
      | (center, members, size) :: centers ->
        begin
          Ht.add clusters center members;
          (* mark members as clustered *)
          L.iter (fun m ->
              Ht.add clustered m ()
            ) members;
          Log.info "cluster: %d size: %d clustered: %d"
            !cid size (Ht.length clustered);
          incr cid;
          let unclustered =
            L.fold_left (fun acc (x, _neighbs, _n) ->
                if Ht.mem clustered x then acc
                else x :: acc
              ) [] centers in
          loop unclustered
        end in
  loop molecules;
  clusters

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
  let threshold = CLI.get_float_def ["-t"] args 0.9 in
  assert(threshold >= 0.0 && threshold <= 1.0);
  verbose := CLI.get_set_bool ["-v"] args;
  CLI.finalize ();
  let all_mols = FpMol.molecules_of_file input_fn in
  let total = L.length all_mols in
  Log.info "read %d from %s" total input_fn;
  Log.info "BST OK";
  let clusters = butina threshold all_mols in
  Log.info "clusters: %d" (Ht.length clusters);
  let i = ref 0 in
  Utls.with_out_file output_fn (fun out ->
      Ht.iter (fun center members ->
          fprintf out "cid: %d size: %d center: %s members:"
            !i (L.length members) (FpMol.get_name center);
          incr i;
          L.iter (fun mol ->
              fprintf out " %s" (FpMol.get_name mol)
            ) members;
          fprintf out "\n"
        ) clusters
    )

let () = main ()
