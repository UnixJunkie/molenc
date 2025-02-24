#!/usr/bin/env python3

# create an HTML table to view all molecules in a web browser
# .smi, .sdf and .mol2 files are supported

import argparse, mols2grid, os, rdkit, sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import rdkit.Chem.rdDepictor
from rdkit.Chem.rdDepictor import ConstrainedDepictionParams

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
    subset = ["img"]
    # <CLI> -------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description = "output a molecules grid to an html file")
    parser.add_argument("-i", metavar = "input_fn", dest = "input_fn",
                        help = "input file {smi|sdf|mol2}")
    parser.add_argument("-o", metavar = "output_fn", dest = "output_fn",
                        help = "output file (html)")
    parser.add_argument("-c", metavar = "num_cols", dest = "n_cols",
                        type = int, default = 5,
                        help = "number of columns (default=5)")
    parser.add_argument('-s', dest='show_names',
                        action='store_true', default=False,
                        help = "show molecule names")
    parser.add_argument('-a', dest='align',
                        action='store_true', default=False,
                        help = "_try_ to depict similarly each pair of molecules")
    # handle options
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    n_cols = args.n_cols
    show_names = args.show_names
    align_2D_drawings = args.align
    # </CLI> ------------------------------------------------------------------
    if show_names:
        subset=["img", "name"]
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

    # try to have pairs of molecules w/ matching 2D drawing/orientation
    # e.g. mol_0 must be aligned to mol_1, mol_2 to mol_3, etc.
    errors = 0
    if align_2D_drawings:
        # precalculate 2D coords for all mols
        for mol in mols:
            AllChem.Compute2DCoords(mol)
        i = 0
        params = ConstrainedDepictionParams()
        params.alignOnly = True
        while i < len(mols) - 1:
            mcs = rdFMCS.FindMCS([mols[i], mols[i+1]])
            mcs_smarts = mcs.smartsString
            mcs_mol = Chem.MolFromSmarts(mcs_smarts)
            # common reference drawing
            AllChem.Compute2DCoords(mcs_mol)
            AllChem.GenerateDepictionMatching2DStructure(mols[i], mcs_mol, params=params)
            AllChem.GenerateDepictionMatching2DStructure(mols[i+1], mcs_mol, params=params)
            i += 2
    if errors > 0:
        print("WARN: superposition errors: %d" % errors, file=sys.stderr)

    # create HTML document
    mols2grid.save(mols,
                   subset = subset,
                   n_cols = n_cols,
                   # size = (300, 350),
                   use_coords = align_2D_drawings,
                   output=output_fn, template="static", prerender=True)

    # view in browser
    cmd = "firefox %s" % output_fn
    print("trying: %s" % cmd, file=sys.stderr)
    os.system(cmd)

if __name__ == '__main__':
    main()
