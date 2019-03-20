
(* Multi-Scale-Encoded molecule *)

module StringMap = BatMap.String

type t = { name: string; map: int StringMap.t }

let create name map =
  { name; map }

let feat_count_of_string s =
  try Scanf.sscanf s "%s %d" (fun s d -> (s, d))
  with exn -> (eprintf "MSE_mol.feat_count_of_string: cannot parse: %s" s;
               raise exn)

let of_lines = function
  | [] -> assert(false)
  | name :: feat_count_strs ->
    let map =
      L.fold_left (fun acc line ->
          let feat, count = feat_count_of_string line in
          StringMap.add feat count acc
        ) StringMap.empty feat_count_strs in
    create name map

let dist x y =
  failwith "not implemented yet"
