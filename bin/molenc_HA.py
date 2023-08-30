#!/usr/bin/env python3

# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

# Robust Heavy Atom count from a .smi input file
# (molenc_lizard.py can ignore some molecules because of RDKit)

import rdkit, sys
from rdkit import Chem
from rdkit.Chem import Lipinski

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for i, line in enumerate(f):
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (Chem.MolFromSmiles(smile, sanitize=False), name)

input_smi = sys.argv[1]
for mol, name in RobustSmilesMolSupplier(input_smi):
    if mol is None:
        print("rdkit could not parse: %s" % name, file=sys.stderr)
    else:
        HA = Lipinski.HeavyAtomCount(mol)
        print("%s\t%d" % (name, HA))
