#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit (CLI tool)

import argparse, rdkit, re, sys, sympy, time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import BRICS, Recap

#TODO
# - remove primes product hack; use a hash table of strings for each
#   atom index instead of atomMapNum; later, those string can be converted
#   back to small numbers

# support up to 100 fragments
frag_identifiers = list(sympy.primerange(2, 541))
num_identifiers = len(frag_identifiers)

digits_star = re.compile('[0-9]+\*')

def remove_BRICS_tags(smi):
    return re.sub(digits_star, '*', smi)

# FBR: return fragments from largest to smallest, to avoid interger overflow
def fragment(recap, mol):
    if recap:
        hierarch = Recap.RecapDecompose(mol)
        return hierarch.GetLeaves().keys()
    else:
        frags = reversed(sorted(BRICS.BRICSDecompose(mol)))
        res = []
        for f in frags:
            res.append(remove_BRICS_tags(f))
        return res

def fragment_RECAP_or_BRICS(recap, out, mol, name):
    frags = fragment(recap, mol)
    num_frags = len(frags)
    scheme = ""
    if recap:
        scheme = "RECAP"
    else:
        scheme = "BRICS"
    print("%d %s fragments in %s" % (num_frags, scheme, name))
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
                            new_id = curr_id * frag_id
                            print(new_id)
                            a.SetAtomMapNum(new_id)
                            a.SetProp("atomNote", str(curr_id * frag_id))
                # a fragment can be matched several times
                # but we want each fragment to have a different index
                frag_idx += 1
    res = Chem.MolToSmiles(mol)
    print('%s\t%s_%s_fragments' % (res, name, scheme), file=out)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules using BRICS or RECAP")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--recap", dest = "recap", action ='store_true',
                        default = True, help = "use RECAP (default)")
    parser.add_argument("--brics", dest = "recap", action ='store_false',
                        help = "use BRICS")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    recap = args.recap
    count = 0
    # fragmentation ---------------------------------------------------------
    with open(output_fn, 'w') as output:
      mol_supplier = RobustSmilesMolSupplier(input_fn)
      for name, mol in mol_supplier:
          print("#atoms:%d %s" % (mol.GetNumAtoms(), name))
          fragment_RECAP_or_BRICS(recap, output, mol, name)
          count += 1
    after = time.time()
    dt = after - before
    print("fragmented %d molecules at %.2f mol/s" % (count, count / dt),
          file=sys.stderr)
