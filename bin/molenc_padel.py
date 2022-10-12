#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Encode molecules using PaDEL molecular descriptors.
#
# Yap, C. W. (2011). PaDEL‐descriptor: An open source software to calculate molecular descriptors and fingerprints. Journal of computational chemistry, 32(7), 1466-1474.

import argparse, padelpy, re, sys

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
    return (smi, name)

def SmilesReader(filename):
    with open(filename) as f:
        for line in f.readlines():
            stripped = line.strip()
            yield parse_smiles_line(stripped)

num_descriptors = 1875 # on 12/10/2022

# line with only default values since PaDEL crashed
def all_missing(mol_name, def_val):
    res = mol_name
    tok = ",%s" % def_val
    for _i in range(1875):
        res += tok
    return res

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Project molecules read from a SMILES file into the \
        PaDEL descriptors space")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file")
    parser.add_argument("-o", metavar = "output_csv", dest = "output_csv",
                        help = "output CSV file")
    parser.add_argument('--no-header', dest='no_header',
                        action='store_true', default=False,
                        help = "no CSV header in output file")
    parser.add_argument('--NA', dest='use_NAs',
                        action='store_true', default=False,
                        help = "use \"NA\" instead of 0 for missing values")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_csv = args.output_csv
    no_header = args.no_header
    use_NAs = args.use_NAs
    default_value = "0"
    if use_NAs:
        default_value = "\"NA\""
    # end CLI parsing ---------------------------------------------------------
    total = 0
    errors = 0
    start = True
    # FBR: we don't handle the case where the 1st molecule would crash PaDEL
    #      in that case we don't even know what the CSV header should be
    with open(output_csv, 'w') as out_file:
        for smi, name in SmilesReader(input_smi):
            descriptors = {}
            try:
                descriptors = padelpy.from_smiles(smi)
            except:
                errors += 1
                print("KO: %s" % smi, file=sys.stderr)
            if (not no_header) and start:
                if errors == 1:
                    print("molenc_padel.py: error on 1st molecule: %s" % smi,
                          file=sys.stderr)
                    sys.exit(1)
                print("\"name\"", end='', file=out_file)
                for k in descriptors:
                    print(",\"%s\"" % k, end='', file=out_file)
                print("", file=out_file) # '\n'
                start = False
            # print all values; start line with molecule name (double quoted)
            to_print = "\"%s\"" % name
            for _k, v in descriptors.items():
                if v == '':
                    # always a joy with molecular descriptors: some
                    # of them cannot be computed for all molecules...
                    to_print += ",%s" % default_value
                else:
                    to_print += ",%g" % float(v)
            if descriptors == {}:
                # robust to PaDEL crashing on some molecules
                to_print = all_missing(name, default_value)
            else:
                assert(len(descriptors) == num_descriptors)
            print(to_print, file=out_file)
            total += 1
    print("encoded: %d errors: %d" % (total, errors), file=sys.stderr)

if __name__ == '__main__':
    main()
