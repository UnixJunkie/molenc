#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# atom typing and molecule fragmentation hints

import argparse, molenc_common, os, random, rdkit, sys, time

from enum import Enum
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Draw import rdMolDraw2D

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
        if ((b.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE) and
            (not b.IsInRing())):
            res.append(b)
    return res

def print_typed_atoms(out, mol):
    for a in mol.GetAtoms():
        i = a.GetIdx()
        t = type_atom(a)
        print("%d %s" % (i, t), file=out)

def char_of_bond_type(bond):
    t = bond.GetBondType()
    if t == rdkit.Chem.rdchem.BondType.SINGLE:
        return '-'
    elif t == rdkit.Chem.rdchem.BondType.AROMATIC:
        return ':'
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
    for bond in bonds:
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        t = char_of_bond_type(bond)
        print("%d %c %d" % (a, t, b), file=out)

# print which bonds are cuttable and the suggested number of cuts
def print_cuttable_bonds(out, mol):
    cuttable_bonds = find_cuttable_bonds(mol)
    total_weight = Descriptors.MolWt(mol)
    # 150 Da: suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cuttable_bonds), nb_frags)
    print("#cut_bonds:%d:%d" % (len(cuttable_bonds), max_cuts), file=out)
    for bond in cuttable_bonds:
        i = bond.GetIdx()
        print("%d" % i, file=out)

# FBR: test on all KEGG drugs
# FBR: the program has two modes: fragment | assemble (TODO)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "compute atom types")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.txt", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--draw", dest = "draw_mol", action ='store_true',
                        default = False,
                        help = "output PNG for each molecule w/ atom indexes")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    draw_mol = args.draw_mol
    output = open(args.output_fn, 'w')
    mol_supplier = RobustSmilesMolSupplier(input_fn)
    count = 0
    for name, mol in mol_supplier:
        print("#atoms:%d %s" % (mol.GetNumAtoms(), name), file=output)
        print_typed_atoms(output, mol)
        print_bonds(output, mol)
        print_cuttable_bonds(output, mol)
        count += 1
        if draw_mol:
            d = rdMolDraw2D.MolDraw2DCairo(500, 500)
            d.drawOptions().addAtomIndices = True
            d.DrawMolecule(mol)
            d.FinishDrawing()
            png_fn = '%s.png' % name
            with open(png_fn, 'wb') as fn:
                fn.write(d.GetDrawingText())
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
