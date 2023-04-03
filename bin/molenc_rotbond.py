#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Compute the diamater of a molecule's 3D conformer
# i.e. largest interatomic distance

import argparse, math, sys
from rdkit import Chem

def euclid(xyz0, xyz1):
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    dx = x0 - x1
    dy = y0 - y1
    dz = z0 - z1
    return math.sqrt(dx*dx + dy*dy + dz*dz)

# WARNING: O(n^2)
def diameter(mol):
    num_atoms = mol.GetNumAtoms()
    conf = mol.GetConformer(0)
    diam = 0.0
    for i in range(num_atoms - 1):
        xyz_i = conf.GetAtomPosition(i)
        for j in range(i + 1, num_atoms):
            xyz_j = conf.GetAtomPosition(j)
            dist = euclid(xyz_i, xyz_j)
            if dist > diam:
                diam = dist
    return diam

if __name__ == '__main__':
    # CLI options parsing
    parser = argparse.ArgumentParser(description =
                                     "compute molecular diameter")
    parser.add_argument("-i", metavar = "input.sdf", dest = "input_fn",
                        help = "3D conformer input file \
                        (single molecule AND conformer)")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    # parse CLI end -----------------------------------------------------------
    count = 0
    mol_supplier = Chem.SDMolSupplier(input_fn)
    for mol in mol_supplier:
        if (mol == None) or (count > 1):
            assert(False)
        count += 1
        diam = diameter(mol)
        print("%f" % diam)
