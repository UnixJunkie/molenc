module Rdkit : sig
  type t

  val of_pyobject : Pytypes.pyobject -> t
  val to_pyobject : t -> Pytypes.pyobject
  val __init__ : smi:string -> unit -> t
  val type_atom : t -> i:int -> unit -> int array
  val get_num_atoms : t -> unit -> int
  val get_distance : t -> i:int -> j:int -> unit -> int
end = struct
  let filter_opt l = List.filter_map Fun.id l

  let py_module =
    lazy
      (let source =
         {pyml_bindgen_string_literal|
import rdkit
from rdkit import Chem

def nb_heavy_atom_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() > 1:
            res += 1
    return res

# return (#HA, #H)
def count_neighbors(a):
    nb_heavy = nb_heavy_atom_neighbors(a)
    nb_H = a.GetTotalNumHs()
    return (nb_heavy, nb_H)

class Rdkit:
    # this is needed because the OCaml side want to know how
    # to get an object of type t
    def __init__(self, smi):
        self.mol = Chem.MolFromSmiles(smi)
        self.mat = Chem.GetDistanceMatrix(self.mol)

    # (atomic_num, #HA, #H, valence - #H, formal_charge)
    def type_atom(self, i):
        a = self.mol.GetAtomWithIdx(i)
        anum = a.GetAtomicNum()
        # do this on the ocaml side, since we get anum
        # assert(anum > 1) # we want to consider only heavy atoms
        nb_HA, nb_H = count_neighbors(a)
        valence = a.GetTotalValence()
        HA_used_val = valence - nb_H
        formal_charge = a.GetFormalCharge()
        return [anum, nb_HA, nb_H, HA_used_val, formal_charge]

    def get_num_atoms(self):
        return self.mol.GetNumAtoms()

    # get the distance (in bonds) between a pair of atoms
    def get_distance(self, i, j):
        return int(self.mat[i][j])
|pyml_bindgen_string_literal}
       in
       let filename =
         {pyml_bindgen_string_literal|rdkit_wrapper.py|pyml_bindgen_string_literal}
       in
       let bytecode = Py.compile ~filename ~source `Exec in
       Py.Import.exec_code_module
         {pyml_bindgen_string_literal|rdkit_wrapper|pyml_bindgen_string_literal}
         bytecode)

  let import_module () = Lazy.force py_module

  type t = Pytypes.pyobject

  let of_pyobject pyo = pyo
  let to_pyobject x = x

  let __init__ ~smi () =
    let callable = Py.Module.get (import_module ()) "Rdkit" in
    let kwargs = filter_opt [ Some ("smi", Py.String.of_string smi) ] in
    of_pyobject @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let type_atom t ~i () =
    let callable = Py.Object.find_attr_string t "type_atom" in
    let kwargs = filter_opt [ Some ("i", Py.Int.of_int i) ] in
    Py.List.to_array_map Py.Int.to_int
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let get_num_atoms t () =
    let callable = Py.Object.find_attr_string t "get_num_atoms" in
    let kwargs = filter_opt [] in
    Py.Int.to_int @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let get_distance t ~i ~j () =
    let callable = Py.Object.find_attr_string t "get_distance" in
    let kwargs =
      filter_opt [ Some ("i", Py.Int.of_int i); Some ("j", Py.Int.of_int j) ]
    in
    Py.Int.to_int @@ Py.Callable.to_function_with_keywords callable [||] kwargs
end
