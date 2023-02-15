open Printf

module A = BatArray
module L = BatList

(* Parallel Gram matrix initialization *)
let emit_one (i: int ref) (n: int) ((): unit): int =
  if !i >= n then raise Parany.End_of_input
  else
    let res = !i in
    incr i;
    res

let process_one (dist: 'a -> 'a -> float) (samples: 'a array) (n: int) (i: int):
  (int * float list) =
  let js = L.range i `To (n - 1) in
  let si = samples.(i) in
  (i, L.map (fun j -> dist si samples.(j)) js)

let gather_one (res: float array array) ((i, xs): (int * float list)): unit =
  L.iteri (fun j' x ->
      let j = j' + i in
      res.(i).(j) <- x;
      res.(j).(i) <- x (* symmetric matrix *)
    ) xs

let initialize_matrix dist ncores csize samples res =
  let n = A.length samples in
  assert(n > 0);
  assert(ncores >= 1);
  if ncores = 1 then (* Sequential *)
    begin
      for i = 0 to n - 1 do
        (* WARNING: we initialize the diagonal while it is all 0s *)
        for j = i to n - 1 do
          let x = dist samples.(i) samples.(j) in
          res.(i).(j) <- x;
          (* WARNING: we could remove the next one *)
          res.(j).(i) <- x (* symmetric matrix *)
        done;
        printf "done: %d/%d\r%!" (i + 1) n;
      done;
      printf "\n%!";
    end
  else (* parallel *)
    Parany.run ~csize ncores
      ~demux:(emit_one (ref 0) n)
      ~work:(process_one dist samples n)
      ~mux:(gather_one res)

(* partial display *)
let print_corners mat =
  let m = A.length mat in
  let n = A.length mat.(0) in
  let idots = ref false in
  for i = 0 to m - 1 do
    if i < 3 || i > m - 4 then
      begin
        let jdots = ref false in
        for j = 0 to n - 1 do
          if j < 3 || j > n - 4 then
            printf (if j <> 0 then "\t%6.2f" else "%6.2f")
              mat.(i).(j)
          else if not !jdots then
            (printf "\t..."; jdots := true)
        done;
        printf "\n"
      end
    else if not !idots then
      (printf "\t\t\t...\n"; idots := true)
  done;
  flush stdout
