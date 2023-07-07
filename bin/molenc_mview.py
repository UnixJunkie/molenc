#!/usr/bin/env python3

# create an HTML table to view all molecules in a web browser
# .smi, .sdf and .mol2 files are supported

import mols2grid, os, rdkit, sys
from rdkit import Chem

input_fn  = sys.argv[1]
output_fn = sys.argv[2]

def sdf_read_mols(fn):
    suppl = Chem.SDMolSupplier(fn)
    return [mol for mol in suppl]

# in a proper .smi file: the first field is the SMILES string
# what's left on the right of it is the molecule's name
def mol_of_smi_line(line):
    strip = line.strip()
    words = strip.split()
    smi = words[0]
    return Chem.MolFromSmiles(smi)

def smi_read_mols(fn):
    lines = open(fn, 'r').readlines()
    return [mol_of_smi_line(line) for line in lines]

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
mols2grid.save(mols, output=output_fn, template="table", prerender=True)

# view in browser
cmd = "firefox %s" % output_fn
print("trying: %s" % cmd, file=sys.stderr)
os.system(cmd)
