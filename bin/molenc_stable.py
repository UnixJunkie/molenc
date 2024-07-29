#!/usr/bin/env python3
#
# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# stable filter: only non-reactive molecules are printed on stdout
# input line format: <SMILES:str>\t<NAME:str>
# output line format: same as input

import sys

from rdkit import Chem
from rdkit.Chem import Descriptors

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            splits = line.strip().split("\t") # enforce TAB-separated
            smile = splits[0]
            try:
                mol = Chem.MolFromSmiles(smile)
                yield (mol, line)
            except Exception:
                print("ERROR: cannot parse: %s" % line,
                      file=sys.stderr)

# Lisurek, M., Rupp, B., Wichard, J., Neuenschwander, M., von Kries, J. P.,
# Frank, R., ... & Kühne, R. (2010).
# Design of chemical libraries with potentially bioactive molecules applying
# a maximum common substructure concept. Molecular diversity, 14(2), 401-408.
# SMARTS patterns kindly provided by Michael Lisurek
pat1 = Chem.MolFromSmarts('[C,c]S(=O)(=O)[F,Cl,Br,I]') # sulfonylhalide
pat2 = Chem.MolFromSmarts('[C,c]S(=O)(=O)O[CX4]') # sulfone_ester
pat3 = Chem.MolFromSmarts('C(=O)[F,Cl,Br,I]') # acylhalide
pat4 = Chem.MolFromSmarts('O=COC=O') # acidanhydride
pat5 = Chem.MolFromSmarts('c1([F,Cl,Br,I])ncccn1') # 2-halo_pyrimidine
pat6 = Chem.MolFromSmarts('[H]C=O') # aldehyde
pat7 = Chem.MolFromSmarts('C(=O)C(=O)') # 1,2-dicarbonyl
pat8 = Chem.MolFromSmarts('C1OC1') # epoxide
pat9 = Chem.MolFromSmarts('C1NC1') # aziridine
pat10 = Chem.MolFromSmarts('C(=O)S') # thioester
pat11 = Chem.MolFromSmarts('[#7]!@[#7]') # hydrazine
pat12 = Chem.MolFromSmarts('C=[CH2]') # ethenes
pat13 = Chem.MolFromSmarts('[H,*,!N][N;!R]=[C;!R]([*,H])[*,H]') # imine
pat14 = Chem.MolFromSmarts('[CX4]I') # alkyl_iodide
pat15 = Chem.MolFromSmarts('[Se]') # selenide
pat16 = Chem.MolFromSmarts('O-O') # peroxide
pat17 = Chem.MolFromSmarts('[NX3]!@[OX2]') # hetero-hetero_single_bond
pat18 = Chem.MolFromSmarts('[NX3]!@[NX3]') # hetero-hetero_single_bond
pat19 = Chem.MolFromSmarts('[NX3]!@[SX2]') # hetero-hetero_single_bond
pat20 = Chem.MolFromSmarts('[SX2]!@[SX2]') # hetero-hetero_single_bond
pat21 = Chem.MolFromSmarts('[SX2]!@[OX2]') # hetero-hetero_single_bond

def stable_filter(mol):
    return (not (
        mol.HasSubstructMatch(pat1) or
        mol.HasSubstructMatch(pat2) or
        mol.HasSubstructMatch(pat3) or
        mol.HasSubstructMatch(pat4) or
        mol.HasSubstructMatch(pat5) or
        mol.HasSubstructMatch(pat6) or
        mol.HasSubstructMatch(pat7) or
        mol.HasSubstructMatch(pat8) or
        mol.HasSubstructMatch(pat9) or
        mol.HasSubstructMatch(pat10) or
        mol.HasSubstructMatch(pat11) or
        mol.HasSubstructMatch(pat12) or
        mol.HasSubstructMatch(pat13) or
        mol.HasSubstructMatch(pat14) or
        mol.HasSubstructMatch(pat15) or
        mol.HasSubstructMatch(pat16) or
        mol.HasSubstructMatch(pat17) or
        mol.HasSubstructMatch(pat18) or
        mol.HasSubstructMatch(pat19) or
        mol.HasSubstructMatch(pat20) or
        mol.HasSubstructMatch(pat21)))

input_fn = sys.argv[1]

for mol, line in RobustSmilesMolSupplier(input_fn):
    if stable_filter(mol):
        # exact input lines replicated to the output
        print(line, end='')
