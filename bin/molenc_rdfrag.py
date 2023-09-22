#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit (CLI tool)

import argparse, rdkit, sys, time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Recap

#TODO
# - find cut bonds; remove them
#   - for each fragment; check it is found only once
# - output fragments as a SMILES mixture

# REMARKS
# FBR: RECAP: is there really a synthesis tree output? With layers?
# FBR: BRICS: how many times a given fragment is matched?

def numerate_atoms(mol):
    i = 1 # atom map nums start at 1
    for a in mol.GetAtoms():
        a.SetAtomMapNum(i)
        i += 1

# # substructure search
# m.HasSubstructMatch(patt)
# True
# m.GetSubstructMatch(patt)
# (0, 5, 6)
# m.GetSubstructMatches(patt)
# ((0, 5, 6), (4, 5, 6))        

def fragment_RECAP(out, mol, name):
    hierarch = Recap.RecapDecompose(mol)
    frags = hierarch.GetLeaves().keys()
    num_frags = len(frags)
    print("fragments: %d" % num_frags)
    if num_frags <= 1:
        print("could not fragment: %s" % Chem.MolToSmiles(mol),
              file=sys.stderr)
    else:
        numerate_atoms(mol)
        in_mol_smi = Chem.MolToSmiles(mol)
        print("%s\t%s" % (in_mol_smi, name), file=out)
        frag_idx = 1
        for smi in frags:
            assert(type(smi) == str)
            patt = Chem.MolFromSmarts(smi)
            matches = mol.GetSubstructMatches(patt)
            num_matches = len(matches)
            if num_matches == 0:
                print("fragment not found: %s" % smi,
                      file=sys.stderr)
            if num_matches > 1:
                print("fragment found several times: %s" % smi,
                      file=sys.stderr)
            for match in matches:
                for i in match:
                    a = mol.GetAtomWithIdx(i)
                    a.SetAtomMapNum(frag_idx)
                frag_idx += 1
    res = Chem.MolToSmiles(mol)
    print('%s\t%s_RECAP_fragments' % (res, name), file=out)

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
