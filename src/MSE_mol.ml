(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* Multi-Scale-Encoded molecule *)

open Printf

module L = MyList
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
