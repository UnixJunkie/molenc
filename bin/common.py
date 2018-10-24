
from __future__ import print_function

import numpy
import rdkit
from rdkit import Chem

# in a bond: atom with lowest index first
# in a list of bonds: bond with lowest first atom index first
def order_bonds_canonically(bonds):
    pairs = map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), bonds)
    min_index_first = map(lambda (a, b): (min(a, b), max(a, b)), pairs)
    min_index_first.sort()
    return min_index_first

def iterate(f, col):
    for x in col:
        f(x)

def print_bond(b):
    print("%d %d" % b)

def print_bonds(mol):
    print("#bonds:%d" % mol.GetNumBonds())
    iterate(print_bond, order_bonds_canonically(mol.GetBonds()))

def print_distance_matrix(mol):
    mat = Chem.GetDistanceMatrix(mol)
    diam = numpy.max(mat)
    print("#diameter:%d" % diam)
    nb_atoms = mol.GetNumAtoms()
    for i in range(nb_atoms):
        for j in range(nb_atoms):
            x = mat[i][j]
            if j == 0:
                print("%d" % x, end='')
            else:
                print(" %d" % x, end='')
        print("") # newline
