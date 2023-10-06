#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit (command-line tool)

import argparse, rdkit, re, sys, sympy, time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import BRICS, Draw, Recap

# set of atom map nums
def atom_map_nums(mol):
    s = set()
    for a in mol.GetAtoms():
        curr_id = a.GetAtomMapNum()
        if curr_id != 0:
            s.add(curr_id)
    return s

# remap elements of an arbitrary integer set w/ cardinality N from 1 to N
def renumber_atom_map_nums(mol):
    s = atom_map_nums(mol)
    # compute the mapping
    ht = {}
    for i, x in enumerate(s):
        ht[x] = i + 1 # AtomMapNums start at 1
    # apply it
    for a in mol.GetAtoms():
        curr_id = a.GetAtomMapNum()
        if curr_id != 0:
            a.SetAtomMapNum(ht[curr_id])

# dump molecule in 2D to a .png file
def png_dump_mol(fn, mol):
    mol_to_draw = Chem.Mol(mol) # copy mol
    renumber_atom_map_nums(mol_to_draw)
    for a in mol_to_draw.GetAtoms():
        curr_id = a.GetAtomMapNum()
        if curr_id != 0:
            # clear AtomMapNum: they look uggly in molecular drawings
            a.SetAtomMapNum(0)
            # use atomNote instead
            a.SetProp("atomNote", str(curr_id))
    print('creating %s' % fn, file=sys.stderr)
    Draw.MolToFile(mol_to_draw, fn, size=(1500,1500))

digits_star = re.compile('[0-9]+\*')

# ignore BRICS "bond tags" for the moment
def remove_BRICS_tags(smi):
    return re.sub(digits_star, '*', smi)

def fragment(recap, mol):
    if recap:
        hierarch = Recap.RecapDecompose(mol)
        return hierarch.GetLeaves().keys()
    else:
        frags = BRICS.BRICSDecompose(mol)
        res = []
        for f in frags:
            res.append(remove_BRICS_tags(f))
        return res

# robust a.GetProp(str_key)
def get_atom_prop(a, key):
    try:
        return a.GetProp(key)
    except KeyError:
        return ""

def dict_value_or_none(d, k):
    try:
        return d[k]
    except KeyError:
        return None

def frag_ids_to_atom_map_nums(mol):
    s = atom_map_nums(mol)
    count = 1
    ht = {}
    # only one of the two fragmenting schemes has already assigned AtomMapNums
    # the other scheme assigns atom properties under the key "frag_ids"
    if len(s) == 0:
        # map atom props to atom map nums
        for a in mol.GetAtoms():
            curr_id = get_atom_prop(a, "frag_ids")
            if curr_id != "":
                val = dict_value_or_none(ht, curr_id)
                if val == None:
                    # new id for this set of fragments
                    ht[curr_id] = count
                    a.SetAtomMapNum(count)
                    count += 1
                else:
                    # already seen set of fragments
                    a.SetAtomMapNum(val)

def fragment_RECAP_or_BRICS(recap, out, mol, name):
    frags = fragment(recap, mol)
    num_frags = len(frags)
    scheme = ""
    if recap:
        scheme = "RECAP"
    else:
        scheme = "BRICS"
    print("%d %s fragments in %s" % (num_frags, scheme, name))
    if num_frags <= 1:
        print("could not fragment: %s" % Chem.MolToSmiles(mol),
              file=sys.stderr)
    else:
        frag_idx = 0
        for smi in frags:
            # /!\ interpreting a SMILES as a SMARTS: is that _really_ safe? /!\
            patt = Chem.MolFromSmarts(smi)
            patt_atoms = patt.GetAtoms()
            matches = mol.GetSubstructMatches(patt)
            if len(matches) == 0:
                print("fragment not found: %s" % smi, file=sys.stderr)
            for match in matches:
                for i, a_idx in enumerate(match):
                    patt_a = patt_atoms[i]
                    # the dummy atom is not part of the fragment
                    if patt_a.GetAtomicNum() > 0:
                        a = mol.GetAtomWithIdx(a_idx)
                        curr_id = get_atom_prop(a, "frag_ids")
                        # don't overwrite already annotated fragments
                        if curr_id == "":
                            a.SetProp("frag_ids", str(frag_idx))
                        else:
                            # atom is part of several fragments
                            # (RECAP is a hierarchical fragmentation scheme)
                            new_id = "%s,%s" % (curr_id, str(frag_idx))
                            a.SetProp("frag_ids", new_id)
                # a fragment can be matched several times
                # but we want each final molecular fragment to have a different id
                frag_idx += 1
    png_fn = '%s.png' % name
    frag_ids_to_atom_map_nums(mol)
    png_dump_mol(png_fn, mol)
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
    print("fragmented %d molecule(s) at %.2f Hz" % (count, count / dt),
          file=sys.stderr)
