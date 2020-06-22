#!/usr/bin/env python3

# One 2D picture SVG for each SMILES line
# molecule images are named after the index of the molecule in the input file
# they are created in the current directory

import argparse
import rdkit
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

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
    for i, name, mol in RobustMolSupplier(input_smi):
        if mol is None:
            continue
        AllChem.Compute2DCoords(mol) # generate 2D conformer
        d = rdMolDraw2D.MolDraw2DSVG(200, 200)
        # d.drawOptions().addAtomIndices = True
        caption = '%d %s' % (i, name)
        d.DrawMolecule(mol, legend = caption)
        d.FinishDrawing()
        out_fn = '%d.svg' % i
        with open(out_fn, 'w') as out:
            out.write(d.GetDrawingText())
