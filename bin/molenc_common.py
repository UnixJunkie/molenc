
import numpy
import rdkit
from rdkit import Chem

def sort_pairs(pairs):
    res = []
    for (a, b) in pairs:
        x = min(a, b)
        y = max(a, b)
        res.append((x,y))
    return res

# in a bond: atom with lowest index first
# in a list of bonds: bond with lowest first atom index first
def order_bonds_canonically(bonds):
    pairs = map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), bonds)
    min_index_first = sort_pairs(pairs)
    min_index_first.sort()
    return min_index_first

def print_bonds(out, mol):
    print("#bonds:%d" % mol.GetNumBonds(), file=out)
    bonds = order_bonds_canonically(mol.GetBonds())
    for b in bonds:
        print("%d %d" % b, file=out)

def print_distance_matrix(out, mol, threeD):
    if threeD:
        # we use a histogram with bin width 1A
        # this allows to work in 3D, at the cost of much more features
        mat = Chem.Get3DDistanceMatrix(mol)
        diam = int(numpy.max(mat))
        print("#diameter:%d" % diam, file=out)
        nb_atoms = mol.GetNumAtoms()
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                x = int(mat[i][j])
                if j == 0:
                    print("%d" % x, end='', file=out)
                else:
                    print(" %d" % x, end='', file=out)
            print("", file=out) # newline
    else:
        mat = Chem.GetDistanceMatrix(mol)
        diam = numpy.max(mat)
        print("#diameter:%d" % diam, file=out)
        nb_atoms = mol.GetNumAtoms()
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                x = mat[i][j]
                if j == 0:
                    print("%d" % x, end='', file=out)
                else:
                    print(" %d" % x, end='', file=out)
            print("", file=out) # newline
