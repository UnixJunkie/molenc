#!/usr/bin/env python

# type atoms of a molecule a la atom pairs
# (nb. pi electrons if > 0, elt. symbol, nbHA neighbors)

from __future__ import print_function

import rdkit
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs

PeriodicTable = Chem.GetPeriodicTable()

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (name, Chem.MolFromSmiles(smile))

def nb_heavy_atom_neighbors(a):
    neighbors = a.GetNeighbors()
    res = 0
    for n in neighbors:
        if n.GetAtomicNum() != 1:
            res = res + 1
    return res

def type_atom(a):
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    symbol = PeriodicTable.GetElementSymbol(a.GetAtomicNum())
    nbHA = nb_heavy_atom_neighbors(a)
    res = ""
    if nb_pi_electrons > 0:
        res = "%d%s%d" % (nb_pi_electrons, symbol, nbHA)
    else:
        res = "%s%d" % (symbol, nbHA)
    return res

def encode_molecule(m):
    atoms = m.GetAtoms()
    res = map(type_atom, atoms)
    return res

def print_encoded_atoms(atoms):
    for a in atoms:
        print(a, end=' ')
    print('\n')

def main():
    if len(sys.argv) != 2:
        print("usage: %s input.smi" % sys.argv[0])
        sys.exit(1)
    input_smi = sys.argv[1]
    # output_csv = sys.argv[2]
    # output = open(output_csv, 'w')
    # output.write(
    #     "#molName logP molMR molW nbA nbD nbRotB TPSA countedAtomPairs...\n")
    for name, mol in RobustSmilesMolSupplier(input_smi):
        if mol is None:
            continue
        print_encoded_atoms(encode_molecule(mol))
    # output.close()

if __name__ == '__main__':
    main()
