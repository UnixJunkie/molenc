(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

   Minimum Spanning Tree (MST) using the dataset's Gram matrix. *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module Utls = Molenc.Utls

(* Undirected graph with integer labels on vertices and float labels (weights)
   on edges *)
module Int = struct
  type t = int
  let compare = BatInt.compare
  let hash = Ht.hash
  let equal = BatInt.equal
  let default = 0
end

module Float = struct
  type t = float
  let compare = BatFloat.compare
  let hash = Ht.hash
  let equal = BatFloat.equal
  let default = 0.0
end

module G = Graph.Imperative.Graph.ConcreteLabeled(Int)(Float)

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

(* FBR: compute the Gram matrix.
 *      Create the graph and populate it with distances from the Gram matrix *)

let mst g =
  Kruskal.spanningtree g

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename>: encode molecules file\n  \
              -o <filename>: output file\n"
       Sys.argv.(0);
     exit 1);
  let _input_fn = CLI.get_string ["-i"] args in
  CLI.finalize ();
  failwith "not implemented yet"

let () = main ()
