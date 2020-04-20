#!/usr/bin/env python3

# InChI to SMILES conversion

import argparse
import rdkit
import sys
from rdkit import Chem

def RobustMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            name = words[0]
            inchi = words[1]
            yield (name, Chem.MolFromInchi(inchi))

if __name__ == '__main__':
    # parse CLI
    # show help in case user has no clue of what to do
    if len(sys.argv) != 3:
        sys.stderr.write("%s input.inchi output.smi\n" % sys.argv[0])
        sys.exit(1)
    input_inchi = sys.argv[1]
    output_smi = sys.argv[2]
    output = open(output_smi, 'w')
    for name, mol in RobustMolSupplier(input_inchi):
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol)
        output.write("%s\t%s\n" % (name, smi))
    output.close()
