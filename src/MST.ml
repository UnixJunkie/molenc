(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   Minimum Spanning Tree (MST) using the dataset's Gram matrix. *)

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
      fprintf out "graph graph_name {\n";
      G.iter_edges (fun src dst ->
          fprintf out "%d -- %d\n" (G.V.label src) (G.V.label dst)
        ) g;
      fprintf out "}\n";
    )

(* write the MST edges to file in dot format *)
let mst_edges_to_dot fn edges =
  Utls.with_out_file fn (fun out ->
      fprintf out "graph graph_name {\n";
      L.iter (fun e ->
          fprintf out "%d -- %d\n" (G.E.src e) (G.E.dst e)
        ) edges;
      fprintf out "}\n";
    )

(* TODO compute the Gram matrix.
 * TODO Create the graph and populate it with distances from the Gram matrix *)

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
    for i = 0 to n - 1 do
      (* WARNING: we initialize the diagonal while it is supposed to be all 0s *)
      for j = i to n - 1 do
        let x = FpMol.dist samples.(i) samples.(j) in
        res.(i).(j) <- x;
        (* WARNING: we could remove the next one *)
        res.(j).(i) <- x (* symmetric matrix *)
      done
    done
  else (* parallel *)
    Parany.run ~csize ncores
      ~demux:(emit_one (ref 0) n)
      ~work:(process_one samples n)
      ~mux:(gather_one res)

(* FBR: output all to all graph to file for inspection/verification *)

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename>: encode molecules file\n  \
              -o <filename>: output file\n  \
              [-np <int>]: maximum number of CPU cores (default=1)\n  \
              [-cs <int>]: parallel job chunk size (default=1)\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
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
  compute_gram_matrix nprocs csize all_mols matrix;
  (* add all edges to graph *)
  for i = 0 to nb_mols - 1 do
    (* WARNING: we don't initialize the diagonal
       (it is supposed to be all 0s) *)
    for j = i + 1 to nb_mols - 1 do
      let w = matrix.(i).(j) in
      let edge = G.E.create i w j in
      G.add_edge_e g edge
    done
  done;
  (* MST *)
  Log.info "MST...";
  let mst = minimum_spanning_tree g in
  (* dump to file *)
  Log.info "dumping to file...";
  mst_edges_to_dot output_fn mst

let () = main ()