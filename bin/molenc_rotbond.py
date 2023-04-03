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

import argparse, sys
from rdkit import Chem

def count_heavy_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() != 1:
            res += 1
    return res

# terminal heavy atom with attached H
def is_hydrogenated_terminal(a):
    return (a.GetAtomicNum() != 1 and         # not H
            count_heavy_neighbors(a) == 1 and # terminal
            a.GetTotalNumHs() >= 1)           # hydrogenated

def already_protonated(mol0):
    before = mol0.GetNumAtoms()
    mol1 = Chem.AddHs(mol0)
    after = mol1.GetNumAtoms()
    return (before == after)

def is_HA(a):
    return (a.GetAtomicNum() != 1)

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
    mol_supplier = Chem.SDMolSupplier(input_fn, removeHs=False)
    for mol in mol_supplier:
        if mol == None:
            assert(False)
        name = mol.GetProp('_Name')
        if not already_protonated(mol):
            print("ERROR: not protonated: %s" % name, file=sys.stderr)
        else:
            regular_bonds = []
            movingH_bonds = []
            for b in mol.GetBonds():
                start_a = b.GetBeginAtom()
                stop_a = b.GetEndAtom()
                # we are only interested in single bonds between heavy atoms
                # and out of rings
                if is_HA(start_a) and is_HA(stop_a) and \
                   b.GetBondTypeAsDouble() == 1.0 and not b.IsInRing():
                    i = start_a.GetIdx()
                    j = stop_a.GetIdx()
                    ij = (i, j)
                    if is_hydrogenated_terminal(start_a) or \
                       is_hydrogenated_terminal(stop_a):
                        movingH_bonds.append(ij)
                    else:
                        regular_bonds.append(ij)
            total = len(regular_bonds) + len(movingH_bonds)
            print('%d:%s' % (total, name))
            for i, j in regular_bonds:
                print("REG\t%d\t%d" % (i, j))
            for i, j in movingH_bonds:
                print("THA\t%d\t%d" % (i, j))
