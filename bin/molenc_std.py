#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# CLI wrapper for rdkit's standardizer

import argparse, re, sys, time
from rdkit import Chem
#from rdkit.Chem.MolStandardize import Standardizer # rdkit.version <= '2023.09.3'
#standardizer = Standardizer()
from rdkit.Chem.MolStandardize import rdMolStandardize # rdkit.version >= '2024.03.5'

regex = re.compile('\\s')

def find_whitespace(s):
    m = re.search(regex, s)
    if m == None:
        return -1
    else:
        return m.start()

def parse_smiles_line(line):
    fst_white = find_whitespace(line)
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
        name = line[fst_white + 1:]
    mol = Chem.MolFromSmiles(smi)
    return (mol, name)

# # for rdkit.version == '2023.09.3'
# def standardize(preserve_stereo, preserve_taut, mol):
#     if preserve_stereo or preserve_taut:
#         s_mol = standardizer.standardize(mol)
#         # We don't need to get fragment parent, because the charge parent is the largest fragment
#         s_mol = standardizer.charge_parent(s_mol, skip_standardize=True)
#         s_mol = standardizer.isotope_parent(s_mol, skip_standardize=True)
#         if not preserve_stereo:
#             s_mol = standardizer.stereo_parent(s_mol, skip_standardize=True)
#         if not preserve_taut:
#             s_mol = standardizer.tautomer_parent(s_mol, skip_standardize=True)
#         return standardizer.standardize(s_mol)
#     else:
#         # standardizer.super_parent(mol): _NOT_ standardizer.standardize(mol)
#         # which doesn't even unsalt the molecule...
#         return standardizer.super_parent(mol)

# for rdkit.version == '2024.03.5'
def standardize(preserve_stereo, preserve_taut, mol):
    if preserve_stereo or preserve_taut:
        # We don't need to get fragment parent, because the charge parent is the largest fragment
        s_mol = rdMolStandardize.ChargeParent(mol, skipStandardize=False)
        s_mol = rdMolStandardize.IsotopeParent(s_mol, skipStandardize=True)
        if not preserve_stereo:
            s_mol = rdMolStandardize.StereoParent(s_mol, skipStandardize=False)
        if not preserve_taut:
            s_mol = rdMolStandardize.TautomerParent(s_mol, skipStandardize=False)
        return s_mol
    else:
        return rdMolStandardize.SuperParent(mol, skipStandardize=False)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "rdkit's standardizer")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output_std.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument('-s', dest='preserve_stereo',
                        action='store_true', default=False,
                        help = "preserve stereo chemistry")
    parser.add_argument('-t', dest='preserve_tautomer',
                        action='store_true', default=False,
                        help = "preserve tautomer (i.e. skip tautomer standardization)")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    preserve_stereo = args.preserve_stereo
    preserve_taut = args.preserve_tautomer
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
                    try:
                        std = standardize(preserve_stereo, preserve_taut, mol)
                        smi_std = Chem.MolToSmiles(std)
                        out.write("%s\t%s\n" % (smi_std, name))
                        count += 1
                    except Exception:
                        print("exception: %s" % name, file=sys.stderr)
                        errors += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; %d errors" % (count, count / dt, errors),
          file=sys.stderr)
