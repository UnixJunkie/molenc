
(* molecular decoder
   output distance of query to each DB mol *)

open Printf

module CLI = Minicli.CLI
module L = MyList
module Ht = BatHashtbl
module StringMap = BatMap.String

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i db\n\
              -i <filename>: encoded molecules database\n"
       Sys.argv.(0);
     exit 1);
  let db_fn = CLI.get_string ["-i"] args in
  let all_lines = Utls.lines_of_file db_fn in
  match all_lines with
  | [] -> assert(false)
  | db_rad :: rest ->
    let all_mols = MSE_mol.of_lines rest in
    printf "%s\n" db_rad;
    let features_seen = Ht.create 11 in
    L.iter (fun mol ->
        let name = MSE_mol.get_name mol in
        let map = MSE_mol.get_map mol in
        printf "%s,0.0,[" name;
        let start = ref true in
        StringMap.iter (fun feat count ->
            let feat_id =
              Ht.(find_default features_seen feat (length features_seen)) in
            Ht.replace features_seen feat feat_id;
            if !start then
              (printf "%d:%d" feat_id count;
               start := false)
            else
              printf ";%d:%d" feat_id count
          ) map;
        printf "]\n"
      ) all_mols;
    Log.info "read %d molecules from %s (%d features)"
      (L.length all_mols) db_fn (Ht.length features_seen)

let () = main ()
