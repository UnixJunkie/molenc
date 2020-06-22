#!/usr/bin/env python3

# One 2D picture SVG for each SMILES line
# molecule images are named after the index of the molecule in the input file
# they are created in the current directory

import argparse
import rdkit
import sys
from rdkit import Chem

def RobustMolSupplier(filename):
    with open(filename) as f:
        i = 0
        for line in f:
            words = line.split()
            index = i
            i += 1
            smi = words[0]
            name = words[1]
            yield (index, name, Chem.MolFromSmiles(smi))

if __name__ == '__main__':
    # parse CLI
    # show help in case user has no clue of what to do
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s input.smi\n" % sys.argv[0])
        sys.exit(1)
    input_smi = sys.argv[1]
    for i, name, mol in RobustMolSupplier(input_inchi):
        if mol is None:
            continue
    output.close()
