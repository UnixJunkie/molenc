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

import argparse
import rdkit
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

from molenc_common import RobustSmilesMolSupplier

def positional_analog_scan(mol, smarts_patt = '[cH]',
                           smi_substs = ['N','CF','CC','CO',
                                         'CCN','CCl','CC(F)(F)(F)','COC']):
    res = []
    ss = set() # a string set
    patt = Chem.MolFromSmarts(smarts_patt)
    for smi in smi_substs:
        subst = Chem.MolFromSmiles(smi)
        analogs = AllChem.ReplaceSubstructs(mol, patt, subst)
        for a in analogs:
            analog_smi = Chem.MolToSmiles(a) # canonicalization
            # remove duplicates
            if analog_smi not in ss:
                res.append(analog_smi)
                ss.add(analog_smi)
    return res

if __name__ == '__main__':
    before = time.time()
    # CLI options
    parser = argparse.ArgumentParser(
        description = "Positional Analog Scanning of each input molecule")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "analogs output file")
    # parse CLI ----------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    count = 0
    # work ----------------------------------------------
    mol_supplier = RobustSmilesMolSupplier(input_fn)
    for name, mol in mol_supplier:
        analogs = positional_analog_scan(mol)
        for i, ana_smi in enumerate(analogs):
            print("%s\t%s_ANA%03d" % (ana_smi, name, i),
                  file=output)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
