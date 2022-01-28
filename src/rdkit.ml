module Rdkit : sig
type t

val of_pyobject : Pytypes.pyobject -> t option

val to_pyobject : t -> Pytypes.pyobject

val __init__ : smi:string -> unit -> t

end = struct
let filter_opt l = List.filter_map Fun.id l

 let import_module () = Py.Import.import_module "rdkit_wrapper" 

type t = Pytypes.pyobject


let is_instance pyo =
  let py_class = Py.Module.get (import_module ()) "Rdkit_wrapper" in
  Py.Object.is_instance pyo py_class

let of_pyobject pyo = if is_instance pyo then Some pyo else None


let to_pyobject x = x

let __init__ ~smi () = let callable = Py.Module.get (import_module ()) "Rdkit_wrapper" in let kwargs = filter_opt [ Some ("smi", Py.String.of_string smi); ] in of_pyobject @@ Py.Callable.to_function_with_keywords callable [||] kwargs

end
