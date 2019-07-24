(* Copyright (C) 2019, Francois Berenger

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
module L = BatList
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
               [-t <float>]: Tanimoto threshold (default=1.0)\n  \
               [-e <filename>]: molecules to exclude from the ones to filter\n  \
               [-a <filename>]: molecules to annotate the ones to filter\n  \
               [-v]: verbose mode\n"
        Sys.argv.(0);
      exit 1
    end;
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let filtered_out_fn = output_fn ^ ".discarded" in
  let threshold = CLI.get_float_def ["-t"] args 1.0 in
  assert(threshold >= 0.0 && threshold <= 1.0);
  verbose := CLI.get_set_bool ["-v"] args;
  (* -a and -e are mutually exclusive options *)
  assert(not (List.mem "-a" args && List.mem "-e" args));
  let mode = match CLI.get_string_opt ["-e"] args with
    | None ->
      begin match CLI.get_string_opt ["-a"] args with
        | None -> Diversify
        | Some fn -> Annotate fn
      end
    | Some fn -> Filter fn in
  CLI.finalize ();
  let threshold_distance = 1.0 -. threshold in
  let read_count = ref 0 in
  let filtered_count = ref 0 in
  match mode with
  | Diversify ->
    let mols_to_filter = FpMol.molecules_of_file input_fn in
    let total = L.length mols_to_filter in
    let ok_mols = diversity_filter threshold_distance mols_to_filter in
    let kept = L.length ok_mols in
    Utls.with_out_file output_fn (fun out ->
        L.iter (fun mol ->
            let name = FpMol.get_name mol in
            fprintf out "%s\n" name
          ) ok_mols
      );
    Log.info "read %d from %s" total input_fn;
    Log.info "kept %d in %s" kept output_fn
  | Filter train_fn ->
    let mols_to_exclude = FpMol.molecules_of_file train_fn in
    let exclude_set = Bstree.of_list mols_to_exclude in
    Utls.with_out_file output_fn (fun out ->
        Utls.with_out_file filtered_out_fn (fun err_out ->
            Utls.iteri_on_lines_of_file input_fn (fun i line ->
                let mol = FpMol.parse_one i line in
                let nearest_train_mol, nearest_d =
                  Bstree.nearest_neighbor mol exclude_set in
                if !verbose then
                  let curr_name = FpMol.get_name mol in
                  let nearest_name = FpMol.get_name nearest_train_mol in
                  printf "%s %s %.2f\n" curr_name nearest_name nearest_d;
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
      )
  | Annotate train_fn ->
    let annot_mols = FpMol.molecules_of_file train_fn in
    let annot_set = Bstree.of_list annot_mols in
    Utls.with_out_file output_fn (fun out ->
        Utls.iteri_on_lines_of_file input_fn (fun i line ->
            let mol = FpMol.parse_one i line in
            let nearest_train_mol, nearest_d =
              Bstree.nearest_neighbor mol annot_set in
            let curr_name = FpMol.get_name mol in
            let nearest_name = FpMol.get_name nearest_train_mol in
            fprintf out "%s %s %.2f\n" curr_name nearest_name nearest_d;
            incr read_count
          )
      );
    Log.info "read %d from %s" !read_count input_fn;
    Log.info "discarded %d in %s" !filtered_count filtered_out_fn

let () = main ()
