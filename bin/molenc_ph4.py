#!/usr/bin/env python3

# project molecules 3D conformers into the pharmacophore features/points space

import os, sys, time
from rdkit import Chem

def matching_indexes(mol, pat_str):
    res = []
    pat = mol.GetSubstructMatches(pat_str)
    for i in pat:
        for j in i:
            res.append(j)
    return res

if __name__ == '__main__':
    before = time.time()
    mol_supplier = Chem.SDMolSupplier(sys.argv[1]) # FBR: handle CLI options properly
    count = 0
    for mol in mol_supplier:
        print("#atoms:%d %s" % (mol.GetNumAtoms(), mol.GetProp('_Name')))
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
