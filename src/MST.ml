(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   Minimum Spanning Tree (MST) over the dataset's Gram matrix.
   Molecules are connected all to all (undirected graph).
   The edge weight between two molecules is the Tanimoto distance between
   their fingerprints.

   Probst, D., & Reymond, J. L. (2020).
   Visualization of very large high-dimensional data sets as minimum spanning
   trees. Journal of Cheminformatics, 12(1), 1-13. *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

(* Undirected graph with integer labels on vertices and float labels (weights)
   on edges *)
module Node = struct
  type t = int
  let compare = BatInt.compare
  let hash = Ht.hash
  let equal = BatInt.equal
  let default = 0
end

module Edge = struct
  type t = float
  let compare = BatFloat.compare
  let hash = Ht.hash
  let equal = BatFloat.equal
  let default = 0.0
end

(* JCF said we should use Graph.Imperative.Matrix.Graph for performance reasons
   because our graphs are dense
   JS said we should use Prim's algorithm instead of Kruskal; also for performance
   reasons since we have a dense graph *)
module G = Graph.Imperative.Graph.ConcreteLabeled(Node)(Edge)

module W = struct
  type label = G.E.label
  type edge = G.E.t
  type t = float (* mandatory *)
  let weight x = G.E.label x
  let zero = 0.0
  let add = (+.)
  let compare = BatFloat.compare (* mandatory *)
end

module Kruskal = Graph.Kruskal.Make(G)(W)

(* write graph to file in graphviz dot format *)
let graph_to_dot fn g =
  Utls.with_out_file fn (fun out ->
      fprintf out "graph all_to_all {\n";
      G.iter_edges_e (fun e ->
          fprintf out "%d -- %d [label=\"%.2f\"]\n"
            (G.E.src e) (G.E.dst e) (G.E.label e)
        ) g;
      fprintf out "}\n";
    )

(* node color *)
let rgb_triplet min_pIC50 delta_pIC50 ic50 =
  let percent = ceil (100.0 *. ((ic50 -. min_pIC50) /. delta_pIC50)) in
  let _fractional, integral = modf percent in
  let key = integral /. 100.0 in
  let red_f, green_f, blue_f = Ht.find Palette.ht key in
  (int_of_float (ceil (255.0 *. red_f)),
   int_of_float (ceil (255.0 *. green_f)),
   int_of_float (ceil (255.0 *. blue_f)))

(* write the MST edges to file in dot format *)
let mst_edges_to_dot fn all_mols edges =
  let pIC50s = A.map FpMol.get_value all_mols in
  let names = A.map FpMol.get_name all_mols in
  let min_pIC50, max_pIC50 = A.min_max pIC50s in
  let delta_pIC50 = max_pIC50 -. min_pIC50 in
  Log.info "pIC50: (min,max,delta): %g %g %g" min_pIC50 max_pIC50 delta_pIC50;
  Utls.with_out_file fn (fun out ->
      fprintf out "graph min_span_tree {\n";
      let nb_nodes = A.length pIC50s in
      for i = 0 to nb_nodes - 1 do
        let ic50 = pIC50s.(i) in
        let name = names.(i) in
        (* color molecule's node *)
        let red, green, blue = rgb_triplet min_pIC50 delta_pIC50 ic50 in
        assert(Utls.in_bounds 0 red   255 &&
               Utls.in_bounds 0 green 255 &&
               Utls.in_bounds 0 blue  255);
        fprintf out "\"%d\" [label=\"\" style=\"filled\" \
                     color=\"#%02x%02x%02x\" \
                     image=\"pix/%s.png\"]\n"
          i red green blue name
      done;
      L.iter (fun e ->
          fprintf out "%d -- %d [label=\"%.2f\"]\n"
            (G.E.src e) (G.E.dst e) (G.E.label e)
        ) edges;
      fprintf out "}\n"
    )

let minimum_spanning_tree g =
  Kruskal.spanningtree g

(* TODO *)
(* FBR: we could use a threshold distance: if two molecules are further than
 *      this distance, we know they are not related (e.g. using DBBAD) *)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename>: encode molecules file\n  \
              -o <filename>: output MST to dot file\n  \
              [-go <filename>]: output fully connected graph to dot file\n  \
              [-np <int>]: maximum number of CPU cores (default=1)\n  \
              [-cs <int>]: parallel job chunk size (default=1)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let maybe_full_graph_fn = CLI.get_string_opt ["-go"] args in
  let nprocs = CLI.get_int_def ["-np"] args 1 in
  let csize = CLI.get_int_def ["-c"] args 1 in
  CLI.finalize ();
  Log.info "reading molecules...";
  let all_mols = A.of_list (FpMol.molecules_of_file input_fn) in
  let nb_mols = A.length all_mols in
  Log.info "read %d" nb_mols;
  let g = G.create ~size:nb_mols () in
  (* add all nodes to graph *)
  for i = 0 to nb_mols - 1 do
    G.add_vertex g i
  done;
  (* compute Gram matrix in // *)
  let matrix = A.init nb_mols (fun _ -> A.create_float nb_mols) in
  Log.info "Gram matrix initialization...";
  Gram.initialize_matrix nprocs csize all_mols matrix;
  Gram.print_corners matrix;
  Log.info "Adding edges to graph...";
  (* add all edges to graph *)
  let disconnected = ref 0 in
  for i = 0 to nb_mols - 1 do
    (* WARNING: we don't initialize the diagonal
       (it is supposed to be all 0s) *)
    for j = i + 1 to nb_mols - 1 do
      let w = matrix.(i).(j) in
      if w = 1.0 then
        incr disconnected
      else
        (* T_dist = 1.0 means the two molecules have nothing in common *)
        (* So, the default molecules "linking threshold" is 1.0. *)
        let edge = G.E.create i w j in
        G.add_edge_e g edge
    done;
    printf "done: %d/%d\r%!" (i + 1) nb_mols;
  done;
  printf "\n%!";
  (if !disconnected > 0 then
     Log.info "disconnected molecules: %d" !disconnected);
  Utls.may_apply (fun fn -> graph_to_dot fn g) maybe_full_graph_fn;
  (* MST *)
  Log.info "MST...";
  let mst = minimum_spanning_tree g in
  (* dump to file *)
  Log.info "writing file %s ..." output_fn;
  mst_edges_to_dot output_fn all_mols mst

let () = main ()
