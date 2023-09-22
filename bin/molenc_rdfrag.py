#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation

import argparse
import molenc_common as common
import rdkit
import sys
import time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "compute molecule fragmentation hints")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.txt", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--draw", dest = "draw_mol", action ='store_true',
                        default = False,
                        help = "output PNG for each molecule w/ atom indexes")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    draw_mol = args.draw_mol
    output = open(args.output_fn, 'w')
    count = 0
    # fragmenting ---------------------------------------------------------
    mol_supplier = RobustSmilesMolSupplier(input_fn)
    for name, mol in mol_supplier:
        print("#atoms:%d %s" % (mol.GetNumAtoms(), name), file=output)
        print_typed_atoms(output, mol)
        print_bonds(output, mol)
        print_cuttable_bonds(output, mol)
        count += 1
        if draw_mol:
            d = rdMolDraw2D.MolDraw2DCairo(500, 500)
            d.drawOptions().addAtomIndices = True
            d.DrawMolecule(mol)
            d.FinishDrawing()
            png_fn = '%s.png' % name
            with open(png_fn, 'wb') as fn:
                fn.write(d.GetDrawingText())
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
