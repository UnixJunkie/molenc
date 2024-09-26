#!/usr/bin/env python3

# create an HTML table to view all molecules in a web browser
# .smi, .sdf and .mol2 files are supported

import argparse, mols2grid, os, rdkit, sys
from rdkit import Chem

def sdf_read_mols(fn):
    suppl = Chem.SDMolSupplier(fn)
    return [mol for mol in suppl]

# in a proper .smi file: the first field is the SMILES string
# what's left on the right of it is the molecule's name
def mol_of_smi_line(line):
    strip = line.strip()
    words = strip.split()
    smi = words[0]
    name = words[1]
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp('name', name)
    return mol

def smi_read_mols(fn):
    lines = open(fn, 'r').readlines()
    return [mol_of_smi_line(line) for line in lines]

def main():
    # <CLI> -------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description = "output molecules grid as html file")
    parser.add_argument("-i", metavar = "input_fn", dest = "input_fn",
                        help = "input file {smi|sdf|mol2}")
    parser.add_argument("-o", metavar = "output_fn", dest = "output_fn",
                        help = "output file (html)")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    # </CLI> ------------------------------------------------------------------
    mols = []
    # ----- SDF -----
    if input_fn.endswith(".sdf"):
        # with some CLI options, we could select a subset:
        # top 50, i..j, last 50, etc.
        mols = sdf_read_mols(input_fn)
    # ----- MOL2 -----
    elif input_fn.endswith(".mol2"):
        print("MOL2: not supported by rdkit; converting to sdf via obabel...", file=sys.stderr)
        input_sdf = input_fn + ".sdf"
        cmd = "obabel %s -O %s" % (input_fn, input_sdf)
        print("trying: %s" % cmd, file=sys.stderr)
        os.system(cmd)
        mols = sdf_read_mols(input_sdf)
    # ----- SMILES -----
    elif input_fn.endswith(".smi"):
        mols = smi_read_mols(input_fn)
    else:
        print("unsupported file type: %s" % input_fn, file=sys.stderr)
        exit(1)

    # create HTML document
    mols2grid.save(mols,
                   # subset=["img", "name"],
                   # n_cols = 2, # designed, nearest_in_training_set
                   # size = (300, 350),
                   output=output_fn, template="static", prerender=True)

    # view in browser
    cmd = "firefox %s" % output_fn
    print("trying: %s" % cmd, file=sys.stderr)
    os.system(cmd)

if __name__ == '__main__':
    main()
