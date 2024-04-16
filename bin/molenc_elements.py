#!/usr/bin/env python3

# output to stdout elements found in each molecule
# INPUT: SMILES file
# OUTPUT: one symbol per line; several times if element present multiple times;
#         hydrogens are made explicit

import sys
from rdkit import Chem

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            smi, _name = line.strip().split("\t") # enforce TAB-separated
            try:
                yield Chem.MolFromSmiles(smi)
            except Exception:
                print("ERROR: cannot parse: %s" % line,
                      file=sys.stderr, end='')

input_fn = sys.argv[1]

for mol in RobustSmilesMolSupplier(input_fn):
    mol_H = Chem.AddHs(mol)
    for a in mol_H.GetAtoms():
        print(a.GetSymbol())
