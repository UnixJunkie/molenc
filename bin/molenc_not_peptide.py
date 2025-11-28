#!/usr/bin/env python3
#
# only output molecules w/o three consecutive peptide bonds
# molecules not filtered out are written to stdout
# usage: molenc_not_peptide.py INPUT.smi

import rdkit
from rdkit import Chem
import sys

input_fn = sys.argv[1]

total = 0
peptides = 0
errors = 0

three_pep_bonds = Chem.MolFromSmiles('NCC(=O)NCC(=O)NCC(=O)')

for line in open(input_fn).readlines():
    smi = line.strip().split()[0]
    mol = Chem.MolFromSmiles(smi)
    if mol:
        if mol.HasSubstructMatch(three_pep_bonds):
            peptides += 1
        else:
            print(line, end='')
    else:
        print('ERROR: could not parse: %s' % line, file=sys.stderr)
        errors += 1
    total += 1

print('INFO: errors/peptides/total=%d/%d/%d' % (errors, peptides, total),
      file=sys.stderr)
