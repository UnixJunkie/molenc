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

(* a blue to white to red color palette *)
(* FBR: this is buggy *)
let rgb_triplet min_pIC50 delta_pIC50 ic50 =
  let fraction = (ic50 -. min_pIC50) /. delta_pIC50 in
  assert(fraction >= 0.0 && fraction <= 1.0);
  if fraction >= 0.5 then
    (* white to red *)
    let red = int_of_float (ceil (255.0 *. fraction)) in
    let green = 255 - red in
    let blue = 255 - red in
    (red, green, blue)
  else (* blue/magenta to white
          (pure blue is too dark to allow reading node labels) *)
    let blue = int_of_float (ceil (255.0 *. (1.0 -. fraction))) in
    let green = blue in
    let red = 255 - green in
    (red, green, blue)

(* write the MST edges to file in dot format *)
let mst_edges_to_dot fn pIC50s edges =
  let min_pIC50, max_pIC50 = A.min_max pIC50s in
  let delta_pIC50 = max_pIC50 -. min_pIC50 in
  Log.info "pIC50: (min,max,delta): %g %g %g" min_pIC50 max_pIC50 delta_pIC50;
  Utls.with_out_file fn (fun out ->
      fprintf out "graph min_span_tree {\n";
      let nb_nodes = A.length pIC50s in
      for i = 0 to nb_nodes - 1 do
        let ic50 = pIC50s.(i) in
        (* color molecule's node *)
        let red, green, blue = rgb_triplet min_pIC50 delta_pIC50 ic50 in
        fprintf out "\"%d\" [label=\"%d\\n%.2f\" \
                     style=\"filled\" color=\"#%02x%02x%02x\"]\n"
          i i ic50 red green blue
      done;
      L.iter (fun e ->
          fprintf out "%d -- %d [label=\"%.2f\"]\n"
            (G.E.src e) (G.E.dst e) (G.E.label e)
        ) edges;
      fprintf out "}\n"
    )

let minimum_spanning_tree g =
  Kruskal.spanningtree g

(* Parallel Gram matrix initialization *)
let emit_one (i: int ref) (n: int) ((): unit): int =
  if !i >= n then raise Parany.End_of_input
  else
    let res = !i in
    incr i;
    res

let process_one (samples: FpMol.t array) (n: int) (i: int):
  (int * float list) =
  let js = L.range i `To (n - 1) in
  let si = samples.(i) in
  (i, L.map (fun j -> FpMol.dist si samples.(j)) js)

let gather_one (res: float array array) ((i, xs): (int * float list)): unit =
  L.iteri (fun j' x ->
      let j = j' + i in
      res.(i).(j) <- x;
      res.(j).(i) <- x (* symmetric matrix *)
    ) xs

let compute_gram_matrix ncores csize samples res =
  let n = A.length samples in
  assert(n > 0);
  assert(ncores >= 1);
  if ncores = 1 then (* Sequential *)
    begin
      for i = 0 to n - 1 do
        (* WARNING: we initialize the diagonal while it is all 0s *)
        for j = i to n - 1 do
          let x = FpMol.dist samples.(i) samples.(j) in
          res.(i).(j) <- x;
          (* WARNING: we could remove the next one *)
          res.(j).(i) <- x (* symmetric matrix *)
        done;
        printf "done: %d/%d\r%!" (i + 1) n;
      done;
      printf "\n%!";
    end
  else (* parallel *)
    Parany.run ~csize ncores
      ~demux:(emit_one (ref 0) n)
      ~work:(process_one samples n)
      ~mux:(gather_one res)

(* FBR: polish the node coloring scheme *)
(* FBR: show the corners of the matrix after init *)
(* FBR: color nodes by percentage relative to (max - min) activity values *)
(* FBR: output the MST in SVG format. Label each node with the SVG of
        the correspondinf 2D molecule *)
(* FBR: we could use a threshold distance: if two molecules are further than
 *      this distance, we know they are not related (e.g. DBBAD) *)

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
  let pIC50s = A.map FpMol.get_value all_mols in
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
  compute_gram_matrix nprocs csize all_mols matrix;
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
  mst_edges_to_dot output_fn pIC50s mst

let () = main ()
