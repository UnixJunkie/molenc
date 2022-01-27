
# type heavy atoms of a molecule and get the distance between each on the molecular graph

import rdkit
from rdkit import Chem

def mol_from_smiles(smi):
    return Chem.MolFromSmiles(smi)

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

# (atomic_num, #HA, #H, valence - #H, formal_charge)
def type_atom(a):
    anum = a.GetAtomicNum()
    assert(anum > 1) # we want to consider only heavy atoms
    nb_HA, nb_H = count_neighbors(a)
    valence = a.GetTotalValence()
    HA_used_val = valence - nb_H
    formal_charge = a.GetFormalCharge()
    return (anum, nb_HA, nb_H, HA_used_val, formal_charge)

# type each atom of the molecule
def type_atoms(mol):
    res = []
    for a in mol.GetAtoms():
        res.append(type_atom(a))
    return res

# get the distance (in bonds) between each pair of atom
def get_distances(mol):
    mat = Chem.GetDistanceMatrix(mol)
    nb_atoms = mol.GetNumAtoms()
    res = []
    for i in range(nb_atoms):
        for j in range(i + 1, nb_atoms):
            dist = int(mat[i][j])
            res.append((i, j, dist))
    return res

# tests
m = Chem.MolFromSmiles('c1ccccc1')
assert(get_distances(m) == [(0, 1, 1),
                            (0, 2, 2),
                            (0, 3, 3),
                            (0, 4, 2),
                            (0, 5, 1),
                            (1, 2, 1),
                            (1, 3, 2),
                            (1, 4, 3),
                            (1, 5, 2),
                            (2, 3, 1),
                            (2, 4, 2),
                            (2, 5, 3),
                            (3, 4, 1),
                            (3, 5, 2),
                            (4, 5, 1)])

assert(type_atoms(m) == [(6, 2, 1, 3, 0),
                         (6, 2, 1, 3, 0),
                         (6, 2, 1, 3, 0),
                         (6, 2, 1, 3, 0),
                         (6, 2, 1, 3, 0),
                         (6, 2, 1, 3, 0)])
