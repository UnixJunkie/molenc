#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Keep only molecules using allowed atoms

import argparse, re, sys, time
from rdkit import Chem

regex = re.compile('\s')

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
        # molecule will percolate instead of behing lost
        smi = line
        name = line
    else:
        smi = line[0:fst_white]
        name = line[fst_white + 1:]
    return Chem.MolFromSmiles(smi)

def parse_atoms_list(line):
    return set(line.strip().split(','))

def atoms_filter(allowed_atoms, mol):
    for a in mol.GetAtoms():
        if a.GetSymbol() not in allowed_atoms:
            return False
    return True

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "filter out molecules w/ disallowed atoms")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument('-a', metavar = "ATOMS_LIST", dest='allowed_atoms',
                        default="C,H,N,O,P,S,F,Cl,Br,I",
                        help = "comma-separated list of allowed atoms \
                        (default=C,H,N,O,P,S,F,Cl,Br,I)")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    allowed_atoms = parse_atoms_list(args.allowed_atoms)
    # parse CLI end -----------------------------------------------------------
    count = 0
    errors = 0
    with open(output_fn, 'w') as out:
        with open(input_fn, 'r') as input:
            for line in input.readlines():
                mol = parse_smiles_line(line.strip())
                if atoms_filter(allowed_atoms, mol):
                    out.write("%s" % line)
                else:
                    errors += 1
                count += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; removed %d" % (count, count / dt, errors),
          file=sys.stderr)
