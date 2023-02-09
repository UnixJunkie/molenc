#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# lead-like filter: only lead-like molecules will be printed on stdout

import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

# Oprea's lead-like filter                                                      
# Hann, M. M., & Oprea, T. I. (2004).                                           
# Pursuing the leadlikeness concept in pharmaceutical research.                 
# Current opinion in chemical biology, 8(3), 255-263.                           
def lead_like(mol):
    # MolW <= 460                                                               
    if Descriptors.MolWt(mol) > 460:
        return False
    # -4.0 <= LogP <= 4.2                                                       
    LogP = Descriptors.MolLogP(mol)
    if LogP < -4.0 or LogP > 4.2:
        return False
    # # LogSw >= -5 # ignored                                                   
    # rotB <= 10                                                                
    if Descriptors.NumRotatableBonds(mol) > 10:
        return False
    # nRings <= 4 (number of SSSR rings, _not_ aromatic rings)                  
    if len(Chem.GetSSSR(mol)) > 4:
        return False
    # HBD <= 5                                                                  
    if Descriptors.NumHDonors(mol) > 5:
        return False
    # HBA <= 9                                                                  
    if Descriptors.NumHAcceptors(mol) > 9:
        return False
    return True # lead-like then!                                               

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
    if lead_like(mol):
        print('%s\t%s' % (smile, name))
