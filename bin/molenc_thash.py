#!/usr/bin/env python3
#
# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# append to stdout <TAB>tautomer_hash to lines of a provided SMILES file

import rdkit, sys, typing
from rdkit import Chem
from rdkit.Chem import RegistrationHash
from rdkit.Chem.RegistrationHash import HashLayer

input_fn = sys.argv[1]

def get_rdkit_tautomer_hash(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    layers = RegistrationHash.GetMolLayers(mol)
    return layers[HashLayer.TAUTOMER_HASH]

for line in open(input_fn).readlines():
    stripped = line.strip()
    smi = stripped.split()[0]
    taut_hash = get_rdkit_tautomer_hash(smi)
    print('%s\t%s' % (stripped, taut_hash))
