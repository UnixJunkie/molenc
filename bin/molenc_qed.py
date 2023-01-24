#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Compute QED for each input SMILES
# scores about 300 molecule/s on a single core
#
# requires: installing qed from sources (pip3 package broken as of 23/01/2023)
# cf. https://github.com/silicos-it/qed
# an open-source implementation of "Quantifying the chemical beauty of drugs"
# https://doi.org/10.1038/nchem.1243

import argparse, rdkit, sys
from qed import qed
from rdkit import Chem

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for i, line in enumerate(f):
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (i, Chem.MolFromSmiles(smile), name)

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Compute Quantitative Estimate of Drug-likeness (QED)")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file")
    parser.add_argument("-o", metavar = "output_tsv", dest = "output_tsv",
                        help = "output CSV file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_tsv = args.output_tsv
    out_count = 0
    error_count = 0
    with open(output_tsv, 'w') as out_file:
        for i, mol, name in RobustSmilesMolSupplier(input_smi):
            if mol is None:
                error_count += 1
            else:
                score = qed.default(mol, False)
                print("%s\t%f" % (name, score), file=out_file)
            out_count += 1
    total_count = out_count + error_count
    print("read: %d errors: %d" % (out_count, error_count),
          file=sys.stderr)

if __name__ == '__main__':
    main()
