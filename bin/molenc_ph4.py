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

# ph4 feature SMARTS from the Pharmer software
# definitions from Lidio Meireles and David Ryan Koes
# Article: https://doi.org/10.1021/ci200097m
# Code: https://raw.githubusercontent.com/UnixJunkie/pharmer/master/pharmarec.cpp

aro_smarts = ["a1aaaaa1",
              "a1aaaa1"]

hbd_smarts = ["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
              "[#8!H0&!$([OH][C,S,P]=O)]",
              "[#16!H0]"]

hba_smarts = ["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
	      "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"]

pos_smarts = ["[+,+2,+3,+4]",
              "[$(CC)](=N)N", # amidine
	      "[$(C(N)(N)=N)]", # guanidine
              "[$(n1cc[nH]c1)]"]

neg_smarts = ["[-,-2,-3,-4]",
              "C(=O)[O-,OH,OX1]",
              "[$([S,P](=O)[O-,OH,OX1])]",
	      "c1[nH1]nnn1",
              "c1nn[nH1]n1",
              "C(=O)N[OH1,O-,OX1]",
              "C(=O)N[OH1,O-]",
	      "CO(=N[OH1,O-])",
              "[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]"] # trifluoromethyl sulfonamide

hyd_smarts = ["a1aaaaa1",
	      "a1aaaa1",
	      # branched terminals as one point
	      "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
	      "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      "*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      # simple rings only; need to combine points to get good results for 3d structures
	      "[C&r3]1~[C&r3]~[C&r3]1",
	      "[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
	      "[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
	      "[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
	      "[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
	      "[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
	      # aliphatic chains
	      "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
	      "[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
	      # sulfur (apparently)
	      "[$([S]~[#6])&!$(S~[!#6])]"]

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
