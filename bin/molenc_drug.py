#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Drug-like filter: only drug-like molecules will be printed on stdout

import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

# Tran-Nguyen, V. K., Jacquemard, C., & Rognan, D. (2020).
# LIT-PCBA: An unbiased data set for machine learning and virtual screening.
# Journal of chemical information and modeling, 60(9), 4263-4273.
def drug_like_filter(mol):
    MolW = Descriptors.MolWt(mol)
    if MolW <= 150 or MolW >= 800: # 150 < MolW < 800 Da
        return False
    cLogP = Descriptors.MolLogP(mol)
    if cLogP <= -3.0 or cLogP >= 5.0: # −3.0 < AlogP < 5.0
        return False
    RotB = Descriptors.NumRotatableBonds(mol)
    if RotB >= 15: # RotB < 15
        return False
    HBA = Descriptors.NumHAcceptors(mol)
    if HBA >= 10: # HBA < 10
        return False
    HBD = Descriptors.NumHDonors(mol)
    if HBD >= 10: # HBD < 10
        return False
    FC = Chem.rdmolops.GetFormalCharge(mol)
    if FC <= -2 or FC >= 2: # −2.0 < FC < 2.0
        return False
    return True # Still here? Drug-like then!

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            smile, name = line.strip().split("\t") # enforce TAB-separated
            try:
                mol = Chem.MolFromSmiles(smile)
                yield (mol, smile, name)
            except Exception:
                print("ERROR: cannot parse: %s" % line,
                      file=sys.stderr, end='')

input_fn = sys.argv[1]

for mol, smile, name in RobustSmilesMolSupplier(input_fn):
    if drug_like_filter(mol):
        print('%s\t%s' % (smile, name))
