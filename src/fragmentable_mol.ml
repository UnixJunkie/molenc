
type atom = { pi_electrons: int;
              atomic_num: int;
              heavy_neighbors: int;
              formal_charge: int }

type bond_type = Single
               | Aromatic
               | Double
               | Triple

type bond = { start: int;
              btype: bond_type;
              stop: int }

(* molecule for which fragmentation hints are known *)
type t = { atoms: atom array;
           bonds: bond array;
           cut_bonds: int array; (* indexes in the bonds array *)
           frag_hint: int } (* how many bonds are suggested to be broken *)
