#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# List rotatable bonds (RotBonds) found in a molecule's 3D conformer.
# Molecules are read from a .sdf file, so they are expected to be in 3D.
# There are two types of rotatable bonds:
# - those that can significantly change the conformation of a molecule
# - those that just make some hydrogens move
#   (hydrogen(s) attached to a "terminal" heavy atom)
# We list them separately, because users might not be interested in trying
# to rotate all rotatable bonds.
# Output format is:
#---
# ^TOTAL_ROT_BONDS:TER_ROT_BONDS:MOL_NAME$
# ^REG	i	j$ # bond_type start_i stop_j
# ...
# ^TER	k	l$
# ...
#---

import argparse, sys
from rdkit import Chem

def already_protonated(mol):
    before = mol.GetNumAtoms()
    molH = Chem.AddHs(mol)
    after = molH.GetNumAtoms()
    return (before == after)

def is_HA(a):
    return (a.GetAtomicNum() != 1)

def count_heavy_neighbors(a):
    res = 0
    for n in a.GetNeighbors():
        if is_HA(n):
            res += 1
    return res

# N with valence 4 becomes N+
def correct_N_formal_charge(mol):
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 7 and \
           a.GetTotalValence() == 4 and \
           a.GetFormalCharge() != +1:
            a.SetFormalCharge(1)
    mol.UpdatePropertyCache()
    Chem.SanitizeMol(mol)
    return mol

# a terminal heavy atom with at least one hydrogen attached
# a.GetTotalNumHs(includeNeighbors=True) avoids an rdkit bug
def is_hydrogenated_terminal(a):
    return (is_HA(a) and \
            count_heavy_neighbors(a) == 1 and \
            a.GetTotalNumHs(includeNeighbors=True) >= 1)

if __name__ == '__main__':
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "list rotatable bonds")
    parser.add_argument("-i", metavar = "input.sdf", dest = "input_fn",
                        help = "3D conformer input file ")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    # parse CLI end -----------------------------------------------------------
    # BUG IF removeHs=False; hydrogens are not properfly counted attached
    # to heavy atoms in this case !!!!!!!!!!!!!!!!
    # REPORTED on ML
    for mol0 in Chem.SDMolSupplier(input_fn, removeHs=False, sanitize=False):
        if mol0 == None:
            assert(False)
        # make rdkit robust to some stupid molecular errors            
        mol = correct_N_formal_charge(mol0)
        name = mol.GetProp('_Name')
        if not already_protonated(mol):
            # a conformer for docking should have explicit hydrogens
            print("ERROR: not protonated: %s" % name, file=sys.stderr)
            # sys.exit(1)
        reg_bonds = []
        ter_bonds = []
        # examine each bond
        for b in mol.GetBonds():
            a_i = b.GetBeginAtom()
            a_j = b.GetEndAtom()
            i = a_i.GetIdx()
            j = a_j.GetIdx()
            ij = (i, j)
            if is_HA(a_i) and is_HA(a_j) and \
               b.GetBondTypeAsDouble() == 1.0 and \
               not b.IsInRing():
                if is_hydrogenated_terminal(a_i) or \
                   is_hydrogenated_terminal(a_j):
                    ter_bonds.append(ij)
                else:
                    reg_bonds.append(ij)
        num_reg_bonds = len(reg_bonds)
        num_ter_bonds = len(ter_bonds)
        total = num_reg_bonds + num_ter_bonds
        print('%d:%d:%s' % (total, num_ter_bonds, name))
        for i, j in reg_bonds:
            print('REG\t%d\t%d' % (i, j))
        for i, j in ter_bonds:
            print('TER\t%d\t%d' % (i, j))
