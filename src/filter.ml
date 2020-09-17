(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Encoded molecules filtering and diversity selection.
   Functionalities:
   - remove training set molecules (or anything too near) from a
     "database" of molecules
     _OR_
   - enforce diversity in a "database" of molecules *)

open Printf

module CLI = Minicli.CLI
module FpMol = Molenc.FpMol
module Ht = Hashtbl
module L = BatList
module Log = Dolog.Log
module TopKeeper = Cpm.TopKeeper
module Utls = Molenc.Utls

module Bstree = struct

  include Bst.Bisec_tree.Make (FpMol)

  let of_list mols =
    create 1 Two_bands (Array.of_list mols)
end

let verbose = ref false

type mode = Filter of string (* file from where to read molecules to exclude;
                                usually, we don't want to order and test
                                molecules which were part of the training set *)
          | Annotate of string (* annotate predicted molecule with nearest
                                  from training set; en easy to understand
                                  explanation for some end users of computer
                                  predictions *)
          | Diversify (* enforce diversity among predicted molecules;
                         might remove some SAR info but OK for lead finding *)
          | K_nearest of int (* retrieve from the "database" the [k]
                                nearest molecules to the query molecule *)

let diversity_filter distance_threshold lst =
  let rec loop acc = function
    | [] -> L.rev acc
    | x :: xs ->
      let name = FpMol.get_name x in
      let ok_mols =
        L.filter (fun y ->
            let d = FpMol.dist x y in
            if d <= distance_threshold then
              let () =
                if !verbose then
                  (* log the ones filtered out and why (the previous molecule
                     which is too near) so that we can inspect them later on *)
                  printf "%s %s %.2f\n" name (FpMol.get_name y) d in
              false
            else
              true
          ) xs in
      loop (x :: acc) ok_mols in
  loop [] lst

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
               [-q <filename>]: file containing a query/reference molecule\n  \
               [-r]: repeat query in output\n  \
               [-k <int>]: number of neighbor molecules to retrieve\n  \
               (required by -q)\n  \
               [-t <float>]: Tanimoto threshold (default=1.0)\n  \
               [-e <filename>]: molecules to exclude from the ones to filter\n  \
               [-a <filename>]: molecules to annotate the ones to filter\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let query_fn = CLI.get_string_opt ["-q"] args in
  let repeat_query = CLI.get_set_bool ["-r"] args in
  let filtered_out_fn = output_fn ^ ".discarded" in
  let threshold = CLI.get_float_def ["-t"] args 1.0 in
  assert(threshold >= 0.0 && threshold <= 1.0);
  verbose := CLI.get_set_bool ["-v"] args;
  (* -a and -e are mutually exclusive options *)
  assert(not (List.mem "-a" args && List.mem "-e" args));
  let mode =
    match query_fn with
    | Some _fn -> (* -q implies -k *)
      let k = CLI.get_int ["-k"] args in
      K_nearest k
    | None ->
      begin match CLI.get_string_opt ["-e"] args with
        | None ->
          begin match CLI.get_string_opt ["-a"] args with
            | None -> Diversify
            | Some fn -> Annotate fn
          end
        | Some fn -> Filter fn
      end in
  CLI.finalize ();
  let threshold_distance = 1.0 -. threshold in
  let read_count = ref 0 in
  let filtered_count = ref 0 in
  match mode with
  | Diversify ->
    begin
      let mols_to_filter = FpMol.molecules_of_file input_fn in
      let total = L.length mols_to_filter in
      let ok_mols = diversity_filter threshold_distance mols_to_filter in
      let kept = L.length ok_mols in
      let ok_names = Ht.create kept in
      L.iter (fun mol ->
              let name = FpMol.get_name mol in
              Ht.add ok_names name ()
            ) ok_mols;
      Utls.with_out_file output_fn (fun out ->
          Utls.iteri_on_lines_of_file input_fn (fun i line ->
              let mol = FpMol.parse_one i line in
              let name = FpMol.get_name mol in
              if Ht.mem ok_names name then
                fprintf out "%s\n" line
            )
        );
      Log.info "read %d from %s" total input_fn;
      Log.info "kept %d in %s" kept output_fn
    end
  | Filter train_fn ->
    begin
      let mols_to_exclude = FpMol.molecules_of_file train_fn in
      let exclude_set = Bstree.of_list mols_to_exclude in
      Utls.with_out_file output_fn (fun out ->
          Utls.with_out_file filtered_out_fn (fun err_out ->
              Utls.iteri_on_lines_of_file input_fn (fun i line ->
                  let mol = FpMol.parse_one i line in
                  let nearest_train_mol, nearest_d =
                    Bstree.nearest_neighbor mol exclude_set in
                  if !verbose then
                    begin
                      let curr_name = FpMol.get_name mol in
                      let nearest_name = FpMol.get_name nearest_train_mol in
                      printf "%s %s %.2f\n" curr_name nearest_name nearest_d
                    end;
                  if nearest_d > threshold_distance then
                    fprintf out "%s\n" line (* keep it *)
                  else
                    begin
                      fprintf err_out "%s\n" line; (* discard it *)
                      incr filtered_count
                    end;
                  incr read_count
                )
            )
        );
      Log.info "read %d from %s" !read_count input_fn;
      Log.info "discarded %d in %s" !filtered_count filtered_out_fn
    end
  | Annotate train_fn ->
    begin
      Log.info "reading %s..." train_fn;
      let annot_mols =
        let mols = FpMol.molecules_of_file train_fn in
        Bstree.of_list mols in
      Log.info "done";
      Utls.with_out_file output_fn (fun out ->
          Utls.iteri_on_lines_of_file input_fn (fun i line ->
              let mol = FpMol.parse_one i line in
              let curr_name = FpMol.get_name mol in
              let nearest_train_mol, nearest_d =
                Bstree.nearest_neighbor mol annot_mols in
              let nearest_name = FpMol.get_name nearest_train_mol in
              let nearest_tani = 1.0 -. nearest_d in
              fprintf out "%s %s %.2f\n" curr_name nearest_name nearest_tani;
              incr read_count;
              Log.info "annotated: %d" !read_count (* user feedback *)
            )
        );
      Log.info "read %d from %s" !read_count input_fn;
      Log.info "discarded %d in %s" !filtered_count filtered_out_fn
    end
  | K_nearest k ->
    begin
      let query_mol = match query_fn with
        | None -> failwith "Filter: K_nearest mode requires -q option"
        | Some fn -> match FpMol.molecules_of_file fn with
          | [q] -> q
          | _ :: _ -> failwith "Filter: several molecules in query file"
          | [] -> failwith "Filter: no molecule in query file"
      in
      let top_k = TopKeeper.create k in
      Utls.with_out_file output_fn (fun out ->
          Utls.iteri_on_lines_of_file input_fn (fun i line ->
              let mol = FpMol.parse_one i line in
              let score = FpMol.tani query_mol mol in
              TopKeeper.add top_k score mol;
              incr read_count
            );
          Log.info "read %d from %s" !read_count input_fn;
          if repeat_query then
            begin
              let query_name = FpMol.get_name query_mol in
              fprintf out "%s_%.3f\n" query_name 1.0
            end;
          let kept = TopKeeper.high_scores_first top_k in
          L.iter (fun (tani, near_mol) ->
              let curr_name = FpMol.get_name near_mol in
              fprintf out "%s_%.3f\n" curr_name tani
            ) kept
        )
    end

let () = main ()
