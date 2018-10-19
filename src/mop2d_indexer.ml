
(* output MOP2D atom environments; sorted by decreasing frequency given
   a database of molecules (a big .mol2 file) *)

open Printf

module L = MyList
module Ht = BatHashtbl
module StringSet = BatSet.String

let main () =
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s -i molecules.mol2 -r max_radius -o output.idx\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let radius = CLI.get_int ["-r"] args in
  assert(BatString.ends_with input_fn ".mol2");
  Utls.with_infile_outfile input_fn output_fn (fun input output ->
      let counter = ref 0 in
      let ht = Ht.create 11 in
      try
        while true do
          let m = Mol2.read_one counter input in
          if !counter mod 1000 = 0 then
            eprintf "%d molecules seen\r%!" !counter; (* user feedback *)
          let mop2d = Mini_mol.mop2d_encode radius m in
          L.iter (fun env ->
              try Ht.modify env ((+) 1) ht
              with Not_found -> Ht.add ht env 1
            ) mop2d
        done
      with End_of_file ->
        let key_values = Ht.to_list ht in
        Log.info "read %d from %s" !counter input_fn;
        let sorted = L.sort (fun (_, v1) (_, v2) ->
            BatInt.compare v2 v1) key_values in
        fprintf output "#radius=%d\n" radius;
        L.iter (fun (k, v) ->
            fprintf output "%s %d\n" (Mop2d_env.to_string k) v
          ) sorted
    )

let () = main ()
