#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# molecule txt fragment to SMILES

import argparse, rdkit, re, sys, time
from rdkit import Chem

# "#atoms:15 NCGC00261552-01_f00"
def read_atoms_header(line):
    (atoms, nb_atoms, frag_name) = [t(s) for t,s in
                                    zip((str,int,str), re.split('[: ]', line))]
    assert(atoms == "#atoms")
    return (nb_atoms, frag_name)

# "0 0,6,2,0"
def read_atom(line):
    (index, nb_pi, atomic_num, nb_HA, charge) = [t(s) for t,s in
                                                 zip((int,int,int,int,int),
                                                     re.split('[, ]', line))]
    return (index, nb_pi, atomic_num, nb_HA, charge)

# "#bonds:16"
def read_bonds_header(line):
    (bonds, nb_bonds) = [t(s) for t,s in
                         zip((str,int), re.split('[:]', line))]
    assert(bonds == "#bonds")
    return nb_bonds

def bond_type_of_char(c):
    if c == '-':
        return rdkit.Chem.rdchem.BondType.SINGLE
    elif c == '~':
        return rdkit.Chem.rdchem.BondType.AROMATIC
    elif c == '=':
        return rdkit.Chem.rdchem.BondType.DOUBLE
    elif c == '#':
        return rdkit.Chem.rdchem.BondType.TRIPLE
    else:
        assert("molenc_frag2smi.py: bond_type_of_char" == "")

# "0 - 16"
def read_bond(line):
    (start_i, c, stop_i) = [t(s) for t,s in
                                          zip((int,str,int),
                                              re.split('[ ]', line))]
    return (start_i, bond_type_of_char(c), stop_i)

# "#anchors:1"
def read_anchors_header(line):
    (anchors, nb_anchors) = [t(s) for t,s in
                             zip((str,int), re.split('[:]', line))]
    assert(anchors == "#anchors")
    return nb_anchors

# "0,6,2,0 0 0,6,2,0"
def read_anchor(line):
    (start_t, start_i, stop_t) = [t(s) for t,s in
                                  zip((str,int,str), re.split('[ ]', line))]
    return start_i

class End_of_file(Exception):
    """End of file was reached"""
    pass

def read_one_fragment(input):
    atoms_header = input.readline()
    if atoms_header == '':
        raise End_of_file
    nb_atoms, frag_name = read_atoms_header(atoms_header)
    atoms = []
    for i in range(nb_atoms):
        line = input.readline()
        atom = read_atom(line)
        atoms.append(atom)
    bonds_header = input.readline()
    nb_bonds = read_bonds_header(bonds_header)
    bonds = []
    for i in range(nb_bonds):
        line = input.readline()
        bond = read_bond(line)
        bonds.append(bond)
    anchors_header = input.readline()
    nb_anchors = read_anchors_header(anchors_header)
    anchors = []
    for i in range(nb_anchors):
        line = input.readline()
        anchor = read_anchor(line)
        anchors.append(anchor)
    print('%s %d %d %d' % (frag_name, nb_atoms, nb_bonds, nb_anchors),
          file=sys.stderr)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "compute atom types")
    parser.add_argument("-i", metavar = "input.frags", dest = "input_fn",
                        help = "fragments input file")
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
                read_one_fragment(input)
                count += 1
        except End_of_file:
            pass
    after = time.time()
    dt = after - before
    print("%d fragments at %.2f frag/s" % (count, count / dt), file=sys.stderr)
    output.close()
