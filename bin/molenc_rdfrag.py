#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit
# a end-user command-line tool

import argparse, rdkit, sys, time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Recap

#TODO
# - find cut bonds; remove them
# - output fragments as a SMILES mixture; keep showing atomMapNums

# REMARKS
# FBR: RECAP: is there really a synthesis tree output? With layers?
# FBR: BRICS: how many times a given fragment is matched?

def numerate_atoms(mol):
    i = 1 # atom map nums start at 1
    for a in mol.GetAtoms():
        a.SetAtomMapNum(i)
        i += 1

def fragment_RECAP(out, mol, name):
    hierarch = Recap.RecapDecompose(mol)
    frags = hierarch.GetLeaves()
    if len(frags) == 1:
        print("could not fragment: %s" % Chem.MolToSmiles(mol),
              file=sys.stderr)
    else:
        numerate_atoms(mol        )
        in_mol_smi = Chem.MolToSmiles(mol)
        print("%s\t%s" % (in_mol_smi, name), file=out)
        # for smi in frags:
        #     print(smi, file=out)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules using BRICS or RECAP")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--brics", dest = "brics", action ='store_true',
                        default = False,
                        help = "use BRICS")
    parser.add_argument("--recap", dest = "recap", action ='store_true',
                        default = True,
                        help = "use RECAP (default)")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    count = 0
    # fragmentation ---------------------------------------------------------
    with open(output_fn, 'w') as output:
      mol_supplier = RobustSmilesMolSupplier(input_fn)
      for name, mol in mol_supplier:
          print("#atoms:%d %s" % (mol.GetNumAtoms(), name))
          fragment_RECAP(output, mol, name)
          count += 1
    after = time.time()
    dt = after - before
    print("fragmented %d molecules at %.2f mol/s" % (count, count / dt),
          file=sys.stderr)
