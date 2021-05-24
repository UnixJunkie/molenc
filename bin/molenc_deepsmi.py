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

from molenc_common import RobustSmilesSupplier

def encode(converter, smi):
    return converter.encode(smi)

def decode(converter, deep_smi):
    try:
        smi = converter.decode(deep_smi)
        # currently, de decoder does not output a canonical SMILES
        # https://github.com/baoilleach/deepsmiles/issues/19
        # I want canonical SMILES, because this is rdkit's default
        mol = Chem.MolFromSmiles(smi)
        cano_smi = Chem.MolToSmiles(mol)
        return cano_smi
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
    parser.add_argument("--no-rings", dest = "rings",
                        action = "store_true",
                        default = False,
                        help = "DeepSMILES without ring openings")
    parser.add_argument("--no-branches", dest = "branches",
                        action = "store_true",
                        default = False,
                        help = "DeepSMILES without branches")
    parser.add_argument("-e", dest = "do_encode",
                        action = "store_true",
                        default = True,
                        help = "encode: SMILES to DeepSMILES (default)")
    parser.add_argument("-d", dest = "do_decode",
                        action = "store_true",
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
    do_encode = args.do_encode
    do_decode = args.do_decode
    if do_decode:
        do_encode = False
    assert(not (do_encode and do_decode))
    if not (rings or branches):
        print("use at least --no-rings or --no-branches",
              file=sys.stderr)
        sys.exit(1)
    count = 0
    # work ----------------------------------------------
    smi_supplier = RobustSmilesSupplier(input_fn)
    converter = deepsmiles.Converter(rings, branches)
    if do_encode:

        for smi, name in smi_supplier:
            deep_smi = encode(converter, smi)
            print("%s\t%s" % (deep_smi, name), file=output)
            count += 1
    else: # decode
        for deep_smi, name in smi_supplier:
            smi = decode(converter, deep_smi)
            if smi != None:
                print("%s\t%s" % (smi, name), file=output)
            count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
