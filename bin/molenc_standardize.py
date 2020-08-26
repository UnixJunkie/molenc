#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.
#
# Molecule standardization and unsalting using rdkit/MolVS

import argparse, rdkit, sys
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import Standardizer

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f.readlines():
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (name, Chem.MolFromSmiles(smile))

salt_remover = SaltRemover()

def unsalt_mol(mol):
    return salt_remover.StripMol(mol, dontRemoveEverything = True)

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Standardize molecules")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file", required = True)
    parser.add_argument("-o", metavar = "output_smi", dest = "output_smi",
                        help = "output SMILES file", required = True)
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_smi = args.output_smi
    error_count = 0
    out_count = 0
    std = Standardizer()
    with open(output_smi, 'w') as out_file:
        for name, mol in RobustSmilesMolSupplier(input_smi):
            if mol is None:
                error_count += 1
            else:
                unsalted = unsalt_mol(mol)
                std_mol = std.standardize(unsalted)
                std_smi = Chem.MolToSmiles(std_mol)
                std_smi_line = "%s\t%s" % (std_smi, name)
                print(std_smi_line, file=out_file)
                out_count += 1
    total_count = out_count + error_count
    print("standardized: %d errors: %d total: %d" % \
          (out_count, error_count, total_count),
          file=sys.stderr)

if __name__ == '__main__':
    main()
