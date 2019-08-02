
open Printf

module CLI = Minicli.CLI
module L = BatList

let main () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let argc, args = CLI.init () in
  if argc = 1 then
    begin
      eprintf "usage:\n\
               %s -i encoded_molecules.txt\n" Sys.argv.(0);
      exit 1
    end;
  let _input_fn = CLI.get_string ["-i"] args in
  let ks = [10; 20; 50; 100] in
  (* read all molecules *)
  L.iter (fun _ k ->
      (* hash them (and compute hashing rate) *)
      (* compute Tani for many pairs (and compute scoring rate) *)
      (* compute estimated tani for the same pairs (and compute scoring rate) *)
      (* output maximum Tani error *)
      failwith "not implemented yet"
    ) ks

let () = main ()
