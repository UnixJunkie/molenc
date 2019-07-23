(* Copyright (C) 2019, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

open Printf

include BatSet.Int

let to_string s =
  let buff = Buffer.create 11 in
  Buffer.add_char buff '[';
  iter (fun x ->
      Buffer.add_string buff
        (if Buffer.length buff = 1
         then sprintf "%d" x
         else sprintf ";%d" x)
    ) s;
  Buffer.add_char buff ']';
  Buffer.contents buff
