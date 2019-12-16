(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Multi-Scale-Encoded molecule *)

open Printf

module L = MyList
module Log = Dolog.Log
module String = BatString
module StringMap = BatMap.String

type t = { name: string; map: int StringMap.t }

let create name map =
  { name; map }

let get_name x =
  x.name

let get_map x =
  x.map

let feat_count_of_string s =
  try Scanf.sscanf s "%s %d" (fun s d -> (s, d))
  with exn -> (eprintf "MSE_mol.feat_count_of_string: cannot parse: %s" s;
               raise exn)

(* to construct one molecules with all its constituent lines
   already read from the input file *)
let read_one = function
  | [] -> failwith "MSE_mol.read_one: empty list"
  | name_line :: feat_count_strs ->
    (* molecule separator is a line starting with a '#' char *)
    assert(String.get name_line 0 = '#');
    let name = String.lchop name_line in (* remove it *)
    let map =
      List.fold_left (fun acc line ->
          let feat, count = feat_count_of_string line in
          (* feature cannot already be here; otherwise,
             there was a problem during encoding of the molecule *)
          if StringMap.mem feat acc then
            Log.warn "mol: %s dup feat: %s" name feat;
          StringMap.add feat count acc
        ) StringMap.empty feat_count_strs in
    create name map

let previous_name = ref ""

exception Break

(* get lines for just one molecule (i.e. for one call to read_one after) *)
let get_lines input =
  let acc = ref [] in
  if !previous_name = "" then
    begin
      let line = input_line input in
      assert(BatString.starts_with line "#"); (* enforce name line *)
      previous_name := line
    end;
  acc := [!previous_name];
  try
    while true do
      let line' = input_line input in
      if BatString.starts_with line' "#" then
        (* this is the start of another molecule *)
        begin
          previous_name := line';
          raise Break
        end
      else
        acc := line' :: !acc
    done;
    assert(false) (* for typing: should never be reached at exec *)
  with Break -> L.rev !acc
     | End_of_file ->
       begin
         previous_name := "";
         L.rev !acc
       end

let of_lines lines =
  let rec loop acc ls =
    match ls with
    | [] -> L.rev acc
    | _ ->
      let name_l, rest =
        L.fold_while (fun l -> String.starts_with l "#")
          (fun acc x -> x :: acc) [] ls in
      (match name_l with
       | [name] ->
         (let feat_counts, remaining_mols =
            L.fold_while (fun l -> not (String.starts_with l "#"))
              (fun acc x -> x :: acc) [] rest in
          let mol = read_one (name :: feat_counts) in
          loop (mol :: acc) remaining_mols)
       | _ -> assert(false)) in
  loop [] lines
