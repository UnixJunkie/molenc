
(* atom type used by the atom pairs molecular encoding; Cf.

   Carhart, R. E., Smith, D. H., & Venkataraghavan, R. (1985).
   Atom pairs as molecular features in structure-activity studies:
   definition and applications.
   Journal of Chemical Information and Computer Sciences, 25(2), 64-73. *)

type t = { pi: int; (* number of pi electrons *)
           elt: string; (* element symbol *)
           nbHA: int } (* number of heavy atom neighbors *)

exception Break

let of_string s =
  let n = String.length s in
  (* extract pi part: leading digits, if any *)
  let nb_leading_digits =
    let count = ref 0 in
    begin
      try
        for i = 0 to n - 1 do
          if BatChar.is_digit s.(i) then incr count
          else raise Break
        done
      with Break -> ()
    end;
    !count in
  (* check format *)
  assert(nb_leading_digits == 0 || nb_leading_digits == 1);
  let pi =
    if nb_leading_digits = 0 then 0
    else Utls.int_of_digit_char s.(0) in
  (* extract elt part: all non digit chars of the string *)
  let elt = BatString.filter (fun c -> not (BatChar.is_digit c)) s in
  (* extract nbHA part: trailing digits, if any *)
  let nb_trailing_digits =
    let count = ref 0 in
    begin
      try
        for i = n - 1 downto 0 do
          if BatChar.is_digit s.(i) then incr count
          else raise Break
        done
      with Break -> ()
    end;
    !count in
  assert(nb_trailing_digits == 0 || nb_trailing_digits == 1);
  let nbHA =
    if nb_trailing_digits = 0 then 0
    else Utls.int_of_digit_char s.(n - 1) in
  { pi; elt; nbHA }
