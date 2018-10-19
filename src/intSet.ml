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
