#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# Only keep molecules with an acceptable number of rotatable bonds

import argparse, re, sys, time
from rdkit import Chem
from rdkit.Chem import Descriptors

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

def rbonds_filter(max_rbonds, mol):
    if Descriptors.NumRotatableBonds(mol) <= max_rbonds:
        return True
    else:
        return False

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "filter out molecules w/ disallowed atoms")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument('-r', metavar = "MAX_ROT_BONDS_INT", dest='max_rbonds',
                        default=-1, type=int,
                        help = "maximum number of rotatable bonds allowed (default=NO_LIMIT")
    # parse CLI ---------------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    max_rbonds = args.max_rbonds
    # parse CLI end -----------------------------------------------------------
    count = 0
    errors = 0
    with open(output_fn, 'w') as out:
        with open(input_fn, 'r') as input:
            for line in input.readlines():
                mol = parse_smiles_line(line.strip())
                if rbonds_filter(max_rbonds, mol):
                    out.write("%s" % line)
                else:
                    errors += 1
                count += 1
    after = time.time()
    dt = after - before
    print("%d molecules @ %.2fHz; removed %d" % (count, count / dt, errors),
          file=sys.stderr)
