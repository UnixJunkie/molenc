#!/usr/bin/env python3

import rdkit, sys
from rdkit import Chem
from rdkit.Chem import Draw

input_smi = sys.argv[1]
output_png = sys.argv[2]

# WARNING: only read and consider first line of input SMILES file
with open(input_smi, 'r') as input:
    with open(output_png, 'wb') as output:
        line = input.readline()
        line.strip()
        split = line.split()
        smi = split[0]
        name = split[1]
        mol = Chem.MolFromSmiles(smi)
        assert(mol != None)
        d2d = Draw.MolDraw2DCairo(-1,-1)
        Draw.DrawMoleculeACS1996(d2d, mol, legend=name)
        pix = d2d.GetDrawingText()
        output.write(pix)
