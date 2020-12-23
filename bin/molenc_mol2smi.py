#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# txt molecule to SMILES

import argparse, rdkit, re, sys, time
import molenc_common as common
from rdkit import Chem

# create a fake molecule for the corresp. fragment
def read_one_molecule(input):
    res_mol = Chem.RWMol()
    atoms_header = input.readline().strip()
    if atoms_header == '':
        raise common.End_of_file # no EOF in Python...
    nb_atoms, name = common.read_atoms_header(atoms_header)
    old2new = {}
    for _i in range(nb_atoms):
        line = input.readline().strip()
        (index, nb_pi, atomic_num, nb_HA, charge) = common.read_atom(line)
        # add atom
        a = Chem.Atom(atomic_num)
        if nb_pi == 1:
            a.SetIsAromatic(True)
        a.SetFormalCharge(charge)
        j = res_mol.AddAtom(a)
        # we need to convert atom indexes
        old2new[index] = j
    bonds_header = input.readline().strip()
    nb_bonds = common.read_bonds_header(bonds_header)
    for i in range(nb_bonds):
        line = input.readline().strip()
        (start_i, bt, stop_i) = common.read_bond(line)
        start = old2new[start_i]
        stop = old2new[stop_i]
        # add bond
        res_mol.AddBond(start, stop, bt)
    # unset aromaticity flag if atom not in ring
    for a in res_mol.GetAtoms():
        if not a.IsInRing():
            a.SetIsAromatic(False)
    try:
        Chem.SanitizeMol(res_mol)
    except rdkit.Chem.rdchem.KekulizeException:
        print("KekulizeException in %s" % name, file=sys.stderr)
    smi = Chem.MolToSmiles(res_mol)
    return (smi, name)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "txt molecule to smi")
    parser.add_argument("-i", metavar = "input.mols", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    count = 0
    with open(input_fn) as input:
        try:
            while True:
                smi, name = read_one_molecule(input)
                count += 1
                print('%s\t%s' % (smi, name), file=output)
        except common.End_of_file:
            pass
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f molecule/s" % (count, count / dt), file=sys.stderr)
    output.close()
