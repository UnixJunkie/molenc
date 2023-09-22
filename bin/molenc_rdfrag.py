#!/usr/bin/env python3
#
# Copyright (C) 2023, Francois Berenger
# BRICS or RECAP molecular fragmentation using rdkit
# (end-user command-line tool)

import argparse
import molenc_common as common
import rdkit
import sys
import time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules using BRICS or RECAP")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--brics", dest = "brics", action ='store_true',
                        default = False,
                        help = "use BRICS")
    parser.add_argument("--recap", dest = "recap", action ='store_true',
                        default = True,
                        help = "use RECAP (default)")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    count = 0
    # fragmenting ---------------------------------------------------------
    with open(output_fn, 'w') as output:
      mol_supplier = RobustSmilesMolSupplier(input_fn)
      for name, mol in mol_supplier:
          print("#atoms:%d %s" % (mol.GetNumAtoms(), name), file=output)
          count += 1
      after = time.time()
      dt = after - before
      print("fragmented %d molecules at %.2f mol/s" % (count, count / dt),
            file=sys.stderr)
