#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Tsuda laboratory, Tokyo University, Japan.

# DeepSMILES encoder/decoder from/to SMILES
#
# "DeepSMILES: An adaptation of SMILES for use in machine-learning of
# chemical structures". Noel M. O’Boyle and Andrew Dalke. ChemRxiv (2018).

import argparse
import deepsmiles
import molenc_common
import rdkit
import sys
import time

from rdkit import Chem
from rdkit.Chem import AllChem

from molenc_common import RobustSmilesMolSupplier

def encode(converter, smi):
    return converter.encode(smi)

def decode(converter, deep_smi):
    try:
        return converter.decode(deep_smi)
    except deepsmiles.DecodeError as e:
        print("molenc_deepsmi.py: decode: '%s'" % e.message,
              file = sys.stderr)
        return None

if __name__ == '__main__':
    before = time.time()
    # CLI options
    parser = argparse.ArgumentParser(
        description = "DeepSMILES encoder/decoder")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument("--no-rings", dest = "rings", action = "store_true",
                        default = False,
                        help = "DeepSMILES without ring openings")
    parser.add_argument("--no-branches", dest = "branches", action = "store_true",
                        default = False,
                        help = "DeepSMILES without branches")
    parser.add_argument("-e", dest = "encode", action = "store_true",
                        default = True,
                        help = "encode: SMILES to DeepSMILES (default)")
    parser.add_argument("-d", dest = "decode", action = "store_true",
                        default = False,
                        help = "decode: DeepSMILES to SMILES")
    # parse CLI ----------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    rings = args.rings
    branches = args.branches
    encode = args.encode
    decode = args.decode
    if encode and decode:
        print("use either -e or -d, not both", file=sys.stderr)
    if (not rings) and (not branches):
        print("use at least --no-rings or --no-branches", file=sys.stderr)
    count = 0
    # work ----------------------------------------------
    converter = deepsmiles.Converter(rings, branches)
    smi_supplier = RobustSmilesSupplier(input_fn)
    if encode:
        for smi, name in smi_supplier:
            deep_smi = encode(converter, smi)
            print("%s\t%s" % (deep_smi, name), file=output)
            count += 1
    else:
        for deep_smi, name in smi_supplier:
            smi = decode(converter, deep_smi)
            if smi != None:
                print("%s\t%s" % (smi, name), file=output)
            count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
