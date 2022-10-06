#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# CLI wrapper for rdkit's standardizer

import argparse, os, string, sys, time
from rdkit import Chem
from rdkit.Chem.MolStandardize import Standardizer

standardizer = Standardizer()

def parse_smiles_line(line):
    fst_space = line.find(' ')
    fst_tab = line.find('\t')
    fst_white = -1
    if fst_space != -1 and fst_tab != -1:
        fst_white = min(fst_space, fst_tab)
    smi = ''
    name = ''
    if fst_white == -1:
        # no whitespace separator: assume molecule has no name
        # use the SMILES itself as the name, so this unnamed
        # molecule will percolate instead instead of behing lost
        smi = line
        name = line
    else:
        smi = line[0:fst_white]
        name = line[fst_white+1:]
    mol = Chem.MolFromSmiles(smi)
    return (mol, name)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "rdkit's standardizer")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output_std.smi", dest = "output_fn",
                        help = "molecules output file")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    # parse CLI end -----------------------------------------------------------
    count = 0
    errors = 0
    with open(output_fn, 'w') as out:
        with open(input_fn, 'r') as input:
            for line in input.readlines():
                stripped = line.strip()
                mol, name = parse_smiles_line(stripped)
                if mol == None:
                    errors += 1
                else:
                    std = s.standardize(mol)
                    smi_std = Chem.MolToSmiles(std)
                    out.write("%s\t%s\n" % (smi_std, name))
                count += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
