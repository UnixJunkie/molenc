(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

open Printf

module A = BatArray
module Ht = BatHashtbl
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

type filename = string

let tap f x =
  f x;
  x

let fst3 (a, _, _) = a

let mkfifo (fn: filename): unit =
  Unix.mkfifo fn 0o600

(* abort if condition is not met *)
let enforce (condition: bool) (err_msg: string): unit =
  if not condition then
    failwith err_msg

let enforce_f (condition: bool) (msg_constr: unit -> string): unit =
  if not condition then
    failwith (msg_constr ())

(* read 'nlines' from 'input' in_channel *)
let read_n_lines nlines input =
  assert(nlines >= 0);
  let rec loop n acc =
    if n = 0 then L.rev acc
    else loop (n - 1) (input_line input :: acc) in
  loop nlines []

(* if the first line is a comment (starts with '#'); then we extract
   it separately from other lines; else we treat it as any other line *)
let maybe_extract_comment_header (fn: filename): string option * string list =
  let all_lines = LO.lines_of_file fn in
  match all_lines with
  | [] -> (None, [])
  | fst :: others ->
    if BatString.starts_with fst "#" then
      (Some fst, others)
    else
      (None, all_lines)

(* Murphy's law: very big files always have problematic lines... *)
let maybe_map_on_lines_of_file (fn: filename) (f: string -> 'a option): 'a list =
  let res = ref [] in
  LO.iter fn (fun line ->
      match f line with
      | Some x -> res := x :: !res
      | None -> Log.warn "line: '%s'" line
    );
  L.rev !res

(* skip 'nb' blocks from file being read *)
let skip_blocks nb read_one input =
  if nb = 0 then ()
  else
    let () = assert(nb > 0) in
    for _ = 1 to nb do
      ignore(read_one input)
    done

(* get the first line (stripped) output by given command *)
let get_command_output (cmd: string): string =
  Log.info "get_command_output: %s" cmd;
  let _stat, output = BatUnix.run_and_read cmd in
  match BatString.split_on_char '\n' output with
  | first_line :: _others -> first_line
  | [] -> (Log.fatal "get_command_output: no output for: %s" cmd; exit 1)

(* run the given command in a sub process (in parallel to the current process)
   and returns its pid so that we can wait for it later *)
let fork_out_cmd (cmd: string): int =
  Log.info "fork_out_cmd: %s" cmd;
  match Unix.fork () with
  | 0 -> (* child process *) exit (Sys.command cmd)
  | -1 -> (* error *) (Log.fatal "fork_out_cmd: fork failed"; exit 1)
  | pid -> pid

(* return full path of command, if found in PATH, none else *)
let command_exists (cmd: string): string option =
  let where_is_cmd = "which " ^ cmd in
  if Unix.system (where_is_cmd ^ " 2>&1 > /dev/null") = Unix.WEXITED 0 then
    Some (get_command_output where_is_cmd)
  else
    None

let run_command (cmd: string): unit =
  Log.info "run_command: %s" cmd;
  match Unix.system cmd with
  | Unix.WSIGNALED _ -> (Log.fatal "run_command: signaled: %s" cmd; exit 1)
  | Unix.WSTOPPED _ -> (Log.fatal "run_command: stopped: %s" cmd; exit 1)
  | Unix.WEXITED i when i <> 0 ->
    (Log.fatal "run_command: exit %d: %s" i cmd; exit 1)
  | Unix.WEXITED _ (* i = 0 then *) -> ()

let get_env (env_var: string): string option =
  try Some (Sys.getenv env_var)
  with Not_found -> None

(* look for exe in PATH then given env. var *)
let find_command (exe: string) (env_var: string): string option =
  match command_exists exe with
  | Some cmd -> Some cmd
  | None ->
    match get_env env_var with
    | Some cmd -> Some cmd
    | None -> (Log.warn "%s not found in PATH; \
                         put it in your PATH or setup the \
                         %s env. var. to point to it" exe env_var;
               None)

let filename_is_absolute fn =
  not (Filename.is_relative fn)

let relative_to_absolute fn =
  if Filename.is_relative fn then
    let cwd = Sys.getcwd () in
    Filename.concat cwd fn
  else
    fn

(* remove the prefix if it is there, or do nothing if it is not *)
let remove_string_prefix prfx str =
  if BatString.starts_with str prfx then
    let prfx_len = String.length prfx in
    BatString.tail str prfx_len
  else
    str

let string_contains_non_binary_digit = Str.regexp "[^01]"

let string_contains_only_zeros_or_ones (s: string): bool =
  not (Str.string_match string_contains_non_binary_digit s 0)

let string_contains_non_digits_non_sep = Str.regexp "[^-0123456789;]"

let string_is_a_list_of_integers (s: string): bool =
  BatString.starts_with s "[" &&
  BatString.ends_with s "]" &&
  let chopped = BatString.chop ~l:1 ~r:1 s in
  not (Str.string_match string_contains_non_digits_non_sep chopped 0)

let may_apply_opt f = function
  | Some x -> Some (f x)
  | None -> None

exception Enough_times

(* accumulate the result of calling 'f' 'n' times *)
let n_times n f =
  let i = ref 0 in
  BatList.unfold_exc (fun () ->
      if !i = n then raise Enough_times
      else
        let res = f () in
        incr i;
        res
    )

(* the identity function *)
let id x = x

let one_or_more_spaces = Str.regexp "[ ]+"

let string_of_floats_array fv =
  let buff = Buffer.create 80 in
  A.iter (fun f ->
      Buffer.add_string buff (sprintf "%f " f)
    ) fv;
  Buffer.contents buff

(* enforce filename uses one of the allowed extensions *)
let enforce_file_extension allowed_exts fn =
  assert(L.exists (BatString.ends_with fn) allowed_exts)

(* Pi math constant *)
let m_pi = 4.0 *. atan 1.0

(* test (lb <= x <= hb) *)
let in_bounds lb x hb =
  x >= lb && x <= hb

let list_medianf (l: float list): float =
  let xs = A.of_list l in
  A.sort BatFloat.compare xs;
  let n = A.length xs in
  if n mod 2 = 1 then
    xs.(n/2)
  else
    0.5 *. (xs.(n/2) +. xs.(n/2 - 1))
(*$T list_medianf
   list_medianf [1.;2.;3.;4.;5.] = 3.0
   list_medianf [1.;2.;3.;4.] = 2.5
*)

let string_of_array ?pre:(pre = "[|") ?sep:(sep = ";") ?suf:(suf = "|]")
    to_str a =
  let buff = Buffer.create 80 in
  Buffer.add_string buff pre;
  A.iteri (fun i x ->
      if i > 0 then Buffer.add_string buff sep;
      Buffer.add_string buff (to_str x)
    ) a;
  Buffer.add_string buff suf;
  Buffer.contents buff

let array_of_string ?pre:(pre = "[|") ?sep:(sep = ";") ?suf:(suf = "|]")
    from_str s =
  let s' = BatString.chop ~l:(String.length pre) ~r:(String.length suf) s in
  let strings = BatString.split_on_string ~by:sep s' in
  let l = L.map from_str strings in
  A.of_list l

(* unmarshal x from file; the file might be gzip, bzip2 or xz compressed *)
let restore fn =
  if S.ends_with fn ".gz" then
    let in_chan = Unix.open_process_in ("gunzip -ck " ^ fn) in
    Marshal.from_channel in_chan
  else if S.ends_with fn ".bz2" then
    let in_chan = Unix.open_process_in ("bunzip2 -ck " ^ fn) in
    Marshal.from_channel in_chan
  else if S.ends_with fn ".xz" then
    let in_chan = Unix.open_process_in ("unxz -ck " ^ fn) in
    Marshal.from_channel in_chan
  else (* assume file is not compressed *)
    LO.restore fn

let is_odd i =
  i mod 2 = 1

let is_even i =
  i mod 2 = 0

let get_first_line fn =
  LO.with_in_file fn input_line

(* like the cut unix command *)
let cut d f line =
  let splitted = BatString.split_on_char d line in
  BatList.at splitted f

let get_ncores () =
  int_of_string (get_command_output "getconf _NPROCESSORS_ONLN")

let int_of_bool = function
  | true -> 1
  | false -> 0

let bool_of_int = function
  | 1 -> true
  | 0 -> false
  | _ -> assert(false)

(* test that 'y' is within 'x' +/- 'epsilon';
   i.e. y \in [x - epsilon, x + epsilon] *)
let approx_equal epsilon x y =
  (y >= x -. epsilon) &&
  (y <= x +. epsilon)

(* proper NaN/nan testing *)
let is_nan x =
  match classify_float x with
  | FP_nan -> true
  | _ -> false

(* some statistics *)

(* population standard deviation *)
let stddev (l: float list): float =
  let n, sx, sx2 =
    List.fold_left (fun (n, sx, sx2) x ->
        (n +. 1., sx +. x, sx2 +. (x *.x))
      ) (0., 0., 0.) l
  in
  sqrt ((sx2 -. (sx *. sx) /. n) /. n)
(* stddev [2.; 4.; 4.; 4.; 5.; 5.; 7.; 9.] = 2.0 *)

(* elements ranks as floats (in a new array) *)
let rank (arr: float array): float array =
  let sorted = A.copy arr in
  A.sort Float.compare sorted;
  let n = A.length arr in
  let ht = Ht.create n in
  let rank = ref 1.0 in
  A.iter (fun x ->
      try
        begin (* deal with ties *)
          let ranks = Ht.find ht x in
          Ht.replace ht x (!rank :: ranks);
          rank := !rank +. 1.0
        end
      with Not_found ->
        begin
          Ht.add ht x [!rank];
          rank := !rank +. 1.0
        end
    ) sorted;
  let elt2rank = Ht.map (fun _elt ranks -> L.favg ranks) ht in
  A.iteri (fun i x ->
      A.unsafe_set sorted i (Ht.find elt2rank x)
    ) arr;
  sorted

(* (\* unit test for rank *\)
 * let () =
 *   assert(rank [|2.; 4.; 7.; 7.; 12.|] = [|1.0; 2.0; 3.5; 3.5; 5.0|]) *)

(* Code comes from Biocaml.Math *)
let wilcoxon_rank_sum_to_z arr1 arr2 =
  let l1, l2 = A.(length arr1, length arr2) in
  let ranked = rank (A.append arr1 arr2) in
  let arr1 = A.sub ranked 0 l1 in
  let l1, l2 = Float.(of_int l1, of_int l2) in
  let sum1 = A.fsum arr1 in
  let expectation = (l1 *. (l1 +. l2 +. 1.)) /. 2. in
  let var = (l1 *. l2 *. ((l1 +. l2 +. 1.) /. 12.)) in
  (sum1 -. expectation) /. (sqrt var)

let cnd x =
  (* Modified from C++ code by David Koppstein. Found from
   www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution *)
  let b1, b2, b3, b4, b5, p, c =
    0.319381530, -0.356563782, 1.781477937, -1.821255978,
    1.330274429, 0.2316419, 0.39894228 in
  if x >= 0.0 then
    let t = 1. /. (1. +. (p *. x)) in
    (1. -. (c *. (exp (-.x *. x /. 2.)) *. t *.
              (t *. (t *. (t *. ((t *. b5) +. b4) +. b3) +. b2) +. b1)))
  else
    let t = 1. /. (1. -. p *. x) in
    c *. (exp (-.x *. x /. 2.)) *. t *.
      (t *. (t *. (t *. ((t *. b5) +. b4) +. b3) +. b2) +. b1)

let wilcoxon_rank_sum_to_p arr1 arr2 =
  (* assumes a two-tailed distribution *)
  let z = wilcoxon_rank_sum_to_z arr1 arr2 in
  2. *. (1. -. (cnd (Float.abs z)))

let fincr_by (xref: float ref) (dx: float): unit =
  xref := !xref +. dx

let int_of_digit_char = function
  | '0' -> 0
  | '1' -> 1
  | '2' -> 2
  | '3' -> 3
  | '4' -> 4
  | '5' -> 5
  | '6' -> 6
  | '7' -> 7
  | '8' -> 8
  | '9' -> 9
  | c -> failwith (sprintf "Utls.int_of_digit_char: not a digit: %c" c)

let char_of_digit = function
  | 0 -> '0'
  | 1 -> '1'
  | 2 -> '2'
  | 3 -> '3'
  | 4 -> '4'
  | 5 -> '5'
  | 6 -> '6'
  | 7 -> '7'
  | 8 -> '8'
  | 9 -> '9'
  | d -> failwith (sprintf "Utls.char_of_digit: not a digit: %d" d)

let string_of_digit d =
  String.make 1 (char_of_digit d)

let make_pair x y =
  (x, y)

(* sort then count duplicates *)
let list_uniq_count l =
  let sorted = L.sort compare l in
  let groups = L.group_consecutive (=) sorted in
  L.map (fun l -> L.(hd l, length l)) groups

(* coarse grain chronometer *)
let time_it f =
  let start = Unix.gettimeofday () in
  let res = f () in
  let stop = Unix.gettimeofday () in
  (stop -. start, res)

let ceili (x: float): int =
  int_of_float (ceil x)

let list_rev_sort cmp l =
  List.sort (fun x y -> cmp y x) l

let array_rand_elt rng a =
  let n = A.length a in
  let i = BatRandom.State.int rng n in
  a.(i)

(* [|1;2;3|] [4;5;6] -> [1;2;3;4;5;6] *)
let prepend_list_with_array a l =
  Array.fold_right (fun x acc ->
      x :: acc
    ) a l

let array_without_elt_at i a =
  let n = (Array.length a) - 1 in
  let res = Array.make n a.(0) in
  let dest = ref 0 in
  for j = 0 to n do
    if j <> i then
      begin
        res.(!dest) <- a.(j);
        incr dest
      end
  done;
  res

let list_really_remove_one l x =
  assert(L.mem x l); (* BatList.remove doesn't enforce this *)
  L.remove l x

let array_prepend_to_list a l =
  A.fold_right (fun x acc -> x :: acc) a l

let string_array_concat sep a =
  let buff = Buffer.create 1024 in
  A.iteri (fun i s ->
      if i > 0 then Buffer.add_string buff sep;
      Buffer.add_string buff s
    ) a;
  Buffer.contents buff

(* approximation of the complementary error function erfc(x),
   comes from the book "numerical recipes" 3rd edition *)
let erfcc x =
  let   z = abs_float x          in
  let   t = 2. /. (2. +. z) in
  let ans = t *. exp
    (-.z *. z -. 1.26551223 +. t *.
              (  1.00002368 +. t *.
              (  0.37409196 +. t *.
              (  0.09678418 +. t *.
              (-.0.18628806 +. t *.
              (  0.27886807 +. t *.
              (-.1.13520398 +. t *.
              (  1.48851587 +. t *.
              (-.0.82215223 +. t *.
                 0.17087277))))))))) in
  if x >= 0.0 then ans
              else 2.0 -. ans

(* Pearson correlation coefficient for float arrays, cross validated with some
   Python implementation of it that I have *)
let pearson_a a1 a2 =
  let tiny = 1.0e-20     in
  let    n = A.length a1 in
  let   nf = float_of_int n      in
  if n <> A.length a2 then
    failwith "my_utils.ml: pearson_a: arrays length differ"
  else
    let p      = A.favg a1 in
    let q      = A.favg a2 in
    let sum_xx = ref 0.       in
    let sum_yy = ref 0.       in
    let sum_xy = ref 0.       in
    let process x' y' =
      let x    = x' -. p in
      let y    = y' -. q in
      let xx   = x *. x  in
      let yy   = y *. y  in
      let xy   = x *. y  in
        sum_xx := !sum_xx +. xx;
        sum_yy := !sum_yy +. yy;
        sum_xy := !sum_xy +. xy;
    in
    for i = 0 to n - 1 do
      process a1.(i) a2.(i);
    done;
    let r = !sum_xy /. (sqrt(!sum_xx *. !sum_yy) +. tiny)        in
    let z = 0.5 *. log((1.0 +. r +. tiny) /. (1.0 -. r +. tiny)) in
    (* approximation of Student's t probability valid for large n *)
    let t = erfcc(abs_float(z *. sqrt(nf -. 1.0)) /. 1.4142136)       in
    (r, t)

(* comes from biocaml, not me *)
let spearman_rank arr =
  let arr = Array.copy arr in
  let arr = Array.mapi (fun i a -> a,i) arr in
  Array.sort (fun (a,_) (b,_) -> BatFloat.compare a b) arr;
  let g _prev il ans =
    let count = List.length il in
    let n = count + (List.length ans) in
    let hi = float_of_int n in
    let lo = float_of_int (n - count + 1) in
    let rank = (hi +. lo) /. 2. in
    (List.map (fun i -> rank,i) il) @ ans
  in
  let f (prev, il, ans) (x,i) =
    let count = List.length il in
    if count = 0 then
      x, [i], ans
    else if x = prev then
      x, i::il, ans
    else
      x, [i], g prev il ans
  in
  let prev,il,ans = Array.fold_left f (0.,[],[]) arr in
  let ans = g prev il ans in
  let ans = List.sort (fun (_,a) (_,b) -> BatInt.compare a b) ans in
  Array.of_list (List.map fst ans)

(* Spearman rank-order correlation coefficient *)
let spearman_a (a1:float array) (a2:float array) =
  pearson_a (spearman_rank a1) (spearman_rank a2)

let spearman_l l1 l2 =
  spearman_a (A.of_list l1) (A.of_list l2)

(* like 'rm -f' *)
let rm_file fn =
  if Sys.file_exists fn then
    Sys.remove fn
