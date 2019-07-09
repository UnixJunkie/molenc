#!/usr/bin/env python3

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor

# FBR: work for a single molecule read from a .smi
#      read atom values from a separate file, one float per line

if __name__ == '__main__':
    m = Chem.MolFromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O") # caffeine
    rdDepictor.Compute2DCoords(m) # create 2D conformer
    selected = [0,1,2,3]
    colors = {0:(1,0,0), 1:(0,1,0), 2:(0,0,1), 3:(1,0,1)}
    drawer = rdMolDraw2D.MolDraw2DSVG(400,400)
    drawer.DrawMolecule(m, highlightAtoms = selected, highlightBonds = [],
                        highlightAtomColors = colors)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    print(svg)
