#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.
#
# Bemis, G. W., & Murcko, M. A. (1996).
# "The properties of known drugs. 1. Molecular frameworks."
# Journal of medicinal chemistry, 39(15), 2887-2893.

import argparse, rdkit, sys
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
        description = "Append Bemis-Murcko scaffold to each input molecule")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file")
    parser.add_argument("-o", metavar = "output_smi", dest = "output_smi",
                        help = "output SMILES file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_smi = args.output_smi
    out_count = 0
    error_count = 0
    with open(output_smi, 'w') as out_file:
        for i, mol, name in RobustSmilesMolSupplier(input_smi):
            if mol is None:
                error_count += 1
            else:
                BM_scaff = 
                print(csv_line, file=out_file)
                out_count += 1
    total_count = out_count + error_count
    print("encoded: %d errors: %d total: %d" % \
          (out_count, error_count, total_count),
          file=sys.stderr)

if __name__ == '__main__':
    main()
