#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# type atoms of a molecule a la atom pairs
# then randomly fragment each input molecule a number of times

import argparse, molenc_common, os, random, rdkit, sys, time

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

# single bonds not in rings
def find_cuttable_bonds(mol):
    res = []
    for b in mol.GetBonds():
        if b.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE and (not b.IsInRing()):
            res.append(b)
    return res

def encode_molecule(m):
    return map(type_atom, m.GetAtoms())

def print_encoded_atoms(out, atoms):
    for i, a in enumerate(atoms):
        print("%d %s" % (i, a), file=out)

def char_of_bond_type(bond):
    print(type(bond))
    t = bond.GetBondType()
    if t == rdkit.Chem.rdchem.BondType.SINGLE:
        return '-'
    elif t == rdkit.Chem.rdchem.BondType.AROMATIC:
        return '~'
    elif t == rdkit.Chem.rdchem.BondType.DOUBLE:
        return '='
    elif t == rdkit.Chem.rdchem.BondType.TRIPLE:
        return '#'
    else:
        assert("molenc_frag.py: char_of_bond_type" == "")

# print all bonds with their type
def print_bonds(out, mol):
    print("#bonds:%d" % mol.GetNumBonds(), file=out)
    bonds = mol.GetBonds()
    for b in bonds:
        a = b.GetBeginAtomIdx()
        b = b.GetEndAtomIdx()
        t = char_of_bond_type(b)
        print("%d %c %d" % (a, t, b), file=out)

# print which bonds are cuttable and the suggested number of cuts
def print_cuttable_bonds(out, mol):
    cuttable_bonds = find_cuttable_bonds(mol)
    total_weight = Descriptors.MolWt(mol)
    # 150Da: suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cuttable_bonds), nb_frags)
    print("#cut_bonds:%d:%d" % (len(cuttable_bonds), max_cuts), file=out)
    for bond in cuttable_bonds:
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        print("%d %d" % (a, b), file=out)

# FBR: test on all KEGG drugs

# FBR: the program has two modes: fragment | assemble

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
        # print_bonds(output, mol)
        print_cuttable_bonds(output, mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
