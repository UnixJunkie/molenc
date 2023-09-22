#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit (CLI tool)

import argparse, rdkit, sys, sympy, time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Recap

#TODO
# - output with more discrete atom numbers (atomMapNum are rendered uggly)

# REMARKS
# FBR: RECAP: is there really a synthesis tree output? With layers?
# FBR: BRICS: how many times a given fragment is matched?

def numerate_atoms(mol):
    i = 1 # atom map nums start at 1
    for a in mol.GetAtoms():
        a.SetAtomMapNum(i)
        i += 1

def clear_atom_map_nums(mol):
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)

# support up to 100 fragments
frag_identifiers = list(sympy.primerange(2, 541))
num_identifiers = len(frag_identifiers)

def fragment_RECAP(out, mol, name):
    hierarch = Recap.RecapDecompose(mol)
    frags = hierarch.GetLeaves().keys()
    num_frags = len(frags)
    print("%d RECAP fragments in %s" % (num_frags, name))
    assert(num_frags < num_identifiers)
    if num_frags <= 1:
        print("could not fragment: %s" % Chem.MolToSmiles(mol),
              file=sys.stderr)
    else:
        # numerate_atoms(mol)
        # in_mol_smi = Chem.MolToSmiles(mol)
        # print("%s\t%s" % (in_mol_smi, name), file=out)
        frag_idx = 0
        for smi in frags:
            #assert(type(smi) == str)
            # print(smi)
            # interpreting a SMILES as a SMARTS: is that _really_ safe ?!
            patt = Chem.MolFromSmarts(smi)
            patt_atoms = patt.GetAtoms()
            matches = mol.GetSubstructMatches(patt)
            if len(matches) == 0:
                print("fragment not found: %s" % smi, file=sys.stderr)
            for match in matches:
                frag_id = frag_identifiers[frag_idx]
                for i, a_idx in enumerate(match):
                    patt_a = patt_atoms[i]
                    # the dummy/* atom is not part of the fragment
                    if patt_a.GetAtomicNum() > 0:
                        a = mol.GetAtomWithIdx(a_idx)
                        curr_id = a.GetAtomMapNum()
                        # don't overwrite already annotated fragments
                        if curr_id == 0:
                            a.SetAtomMapNum(frag_id)
                        else:
                            # atom is part of several fragments
                            # (RECAP is a hierarchical fragmentation scheme)
                            a.SetAtomMapNum(curr_id * frag_id)
                            a.SetProp("atomNote", str(curr_id * frag_id))
                # a fragment can be matched several times
                # but we want each fragment to have a different index
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
