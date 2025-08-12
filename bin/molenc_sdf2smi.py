#!/usr/bin/env python3

import argparse

from rdkit import Chem
from rdkit.Chem.rdmolfiles import SmilesWriter

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help="sdf input file")
parser.add_argument('outputfile', help="smi output file")
args = parser.parse_args()
sdf = Chem.SDMolSupplier(args.inputfile)
writer = SmilesWriter(args.outputfile, delimiter='\t', includeHeader=False)

for mol in sdf:
  if mol is not None:
    writer.write(mol)
writer.close()
