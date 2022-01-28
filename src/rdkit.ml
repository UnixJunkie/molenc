module Rdkit : sig
type t

val of_pyobject : Pytypes.pyobject -> t

val to_pyobject : t -> Pytypes.pyobject

val __init__: smi:string -> unit -> t

val type_atom: t -> i:int -> unit -> int array

val get_num_atoms: t -> unit -> int

val get_distance: t -> i:int -> j:int -> unit -> int

end = struct
let filter_opt l = List.filter_map Fun.id l

 let import_module () = Py.Import.import_module "rdkit_wrapper" 

type t = Pytypes.pyobject

let of_pyobject pyo = pyo

let to_pyobject x = x

let __init__ ~smi () = let callable = Py.Module.get (import_module ()) "Rdkit_wrapper" in let kwargs = filter_opt [ Some ("smi", Py.String.of_string smi); ] in of_pyobject @@ Py.Callable.to_function_with_keywords callable [||] kwargs

let type_atom t ~i () = let callable = Py.Object.find_attr_string t "type_atom" in let kwargs = filter_opt [ Some ("i", Py.Int.of_int i); ] in Py.List.to_array_map Py.Int.to_int @@ Py.Callable.to_function_with_keywords callable [||] kwargs

let get_num_atoms t () = let callable = Py.Object.find_attr_string t "get_num_atoms" in let kwargs = filter_opt [ ] in Py.Int.to_int @@ Py.Callable.to_function_with_keywords callable [||] kwargs

let get_distance t ~i ~j () = let callable = Py.Object.find_attr_string t "get_distance" in let kwargs = filter_opt [ Some ("i", Py.Int.of_int i); Some ("j", Py.Int.of_int j); ] in Py.Int.to_int @@ Py.Callable.to_function_with_keywords callable [||] kwargs

end
