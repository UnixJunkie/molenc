
type t = Single of int (* encoding radius *)
       | Multi of int * int (* start-stop encoding radii *)

let of_string s =
  if BatString.contains s ':' then
    let istr, jstr = BatString.split s ~by:":" in
    Multi (int_of_string istr, int_of_string jstr)
  else
    Single (int_of_string s)
