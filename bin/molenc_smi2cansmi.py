#!/usr/bin/env python3

# SMILES to RdKit canonical SMILES

import sys

from rdkit import Chem

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            smile, name = line.strip().split("\t") # enforce TAB-separated
            try:
                mol = Chem.MolFromSmiles(smile)
                cano_smi = Chem.MolToSmiles(mol)
                yield (cano_smi, name)
            except Exception:
                print("ERROR: cannot parse: %s" % line,
                      file=sys.stderr, end='')

input_fn = sys.argv[1]

for cano_smi, name in RobustSmilesMolSupplier(input_fn):
    print('%s\t%s' % (cano_smi, name))
