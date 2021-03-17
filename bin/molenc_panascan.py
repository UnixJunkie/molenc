#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# Implementation of
# "Positional Analogue Scanning: An Effective Strategy for
# Multiparameter Optimization in Drug Design".
# Pennington, L. D., Aquila, B. M., Choi, Y., Valiulin, R. A., & Muegge, I.
# Journal of medicinal chemistry (2020).
# https://doi.org/10.1021/acs.jmedchem.9b02092

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def positional_analog_scan(mol, smarts_patt = '[cH]',
                           smi_substs = ['N','CF','CC','CO',
                                         'CCN','CCl','CC(F)(F)(F)','COC']):
    res = []
    patt = Chem.MolFromSmarts(smarts_patt)
    for smi in smi_substs:
        subst = Chem.MolFromSmiles(smi)
        analogs = AllChem.ReplaceSubstructs(mol, patt, subst)
        for a in analogs:
            analog_smi = Chem.MolToSmiles(a) # canonicalization
            if not analog_smi in res: # remove duplicates
                res.append(analog_smi)
    return res

# FBR: implement -i and -o
# rename molecules: append _AN%03d to the name

# test
chinine_smi = \
'[H][C@@]1([C@@H](C2=CC=NC3=CC=C(C=C23)OC)O)C[C@@H]4CC[N@]1C[C@@H]4C=C'
print('%s\tchinine' % chinine_smi)
mol = Chem.MolFromSmiles(chinine_smi)
analogs = positional_analog_scan(mol)
for i, a in enumerate(analogs):
    print('%s\tchinine_%02d' % (a, i))
