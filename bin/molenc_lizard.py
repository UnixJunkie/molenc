#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.
#
# Simple molecular encoder using a handful of molecular descriptors.
# type signature: molecule -> (MolW,cLogP,TPSA,RotB,HBA,HBD,FC)
# Why is this script called molenc_lizard.py
# First, because it is part of the molenc project.
# Second, because if you have a little bit of imagination, a lizard
# is just a very small dragon.
#
# The --remove-aliens option was inspired by
# Tran-Nguyen, V. K., Jacquemard, C., & Rognan, D. (2020).
# "LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening."
# Journal of Chemical Information and Modeling.
# https://doi.org/10.1021/acs.jcim.0c00155

#from __future__ import print_function # FBR: relic of python2 ???

import argparse, rdkit, sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs

PeriodicTable = Chem.GetPeriodicTable()

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (Chem.MolFromSmiles(smile), name)

# FBR: add --remove-aliens CLI option

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Project molecules read from a SMILES file into a 7D \
        space whose dimensions are molecular descriptors: \
        (MolW, cLogP, TPSA, RotB, HBA, HBD, FC)")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi")
    parser.add_argument("-o", metavar = "output_csv", dest = "output_csv")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    print("#name\tlogP\tMR\tMW\tHBA\HBD\RotB\tTPSA")
    with open(output_csv, 'w') as out_file:
        for mol, name in RobustSmilesMolSupplier(input_smi):
            if mol is not None:
                molW = Descriptors.MolWt(mol)
                cLogP = Descriptors.MolLogP(mol)
                molMR = Descriptors.MolMR(mol)
                tpsa = Descriptors.TPSA(mol)
                nbRotB = Descriptors.NumRotatableBonds(mol)
                nbA = Descriptors.NumHAcceptors(mol)
                nbD = Descriptors.NumHDonors(mol)
                fc = Chem.rdmolops.GetFormalCharge(mol)
                print("%s %f %f %f %d %d %d %f" %
                      (name, logP, molMR, molW, nbA, nbD, nbRotB, tpsa),
                      file=output_csv)

if __name__ == '__main__':
    main()
