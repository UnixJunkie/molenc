#!/usr/bin/env python3

# One 2D picture SVG for each SMILES line
# molecule images are created in a created pix/ directory
# and named after their corresponding molecule

import argparse
import rdkit
import os
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
    if not (os.path.isdir('pix')):
        os.mkdir('pix')
    for i, name, mol in RobustMolSupplier(input_smi):
        if mol is None:
            continue
        AllChem.Compute2DCoords(mol) # generate 2D conformer
        d = rdMolDraw2D.MolDraw2DCairo(300, 300) # PNG output
        # d.drawOptions().addAtomIndices = True
        caption = '%d %s' % (i, name)
        d.DrawMolecule(mol, legend = caption)
        d.FinishDrawing()
        out_fn = 'pix/%s.png' % name
        print("creating %s" % out_fn)
        with open(out_fn, 'wb') as out:
            out.write(d.GetDrawingText())
