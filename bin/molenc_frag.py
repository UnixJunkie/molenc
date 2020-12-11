#!/usr/bin/env python3

# type atoms of a molecule a la atom pairs
# then fragment each input molecule a number of times

import argparse, molenc_common, os, rdkit, sys, time
from enum import Enum
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = " ".join(words[1:]) # everything after the SMILES string
            yield (name, Chem.MolFromSmiles(smile))

def nb_heavy_atom_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() != 1:
            res += 1
    return res

def type_atom(a):
    # stereo chemistry is ignored for the moment
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    atom_num = a.GetAtomicNum()
    nbHA = nb_heavy_atom_neighbors(a)
    formal_charge = a.GetFormalCharge()
    # make this easy to parse / unambiguous
    res = "%d,%d,%d,%d" % (nb_pi_electrons, atom_num, nbHA, formal_charge)
    return res

def encode_molecule(m):
    return map(type_atom, m.GetAtoms())

def print_encoded_atoms(out, atoms):
    for i, a in enumerate(atoms):
        print("%d %s" % (i, a), file=out)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "compute atom types")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.txt", dest = "output_fn",
                        help = "output file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    mol_supplier = RobustSmilesMolSupplier(input_fn)
    count = 0
    for name, mol in mol_supplier:
        print("#atoms:%d %s" % (mol.GetNumAtoms(), name), file=output)
        print_encoded_atoms(output, encode_molecule(mol))
        molenc_common.print_bonds(output, mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
