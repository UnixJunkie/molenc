#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# atom typing and molecule fragmentation hints

import argparse, molenc_common, os, random, rdkit, re, sys, time

from enum import Enum
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Draw import rdMolDraw2D

# def type_atom(a):
#     # stereo chemistry is ignored for the moment
#     nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
#     atom_num = a.GetAtomicNum()
#     nbHA = nb_heavy_atom_neighbors(a)
#     formal_charge = a.GetFormalCharge()
#     # make this easy to parse / unambiguous
#     res = "%d,%d,%d,%d" % (nb_pi_electrons, atom_num, nbHA, formal_charge)
#     return res

# def char_of_bond_type(bond):
#     t = bond.GetBondType()
#     if t == rdkit.Chem.rdchem.BondType.SINGLE:
#         return '-'
#     elif t == rdkit.Chem.rdchem.BondType.AROMATIC:
#         return '~'
#     elif t == rdkit.Chem.rdchem.BondType.DOUBLE:
#         return '='
#     elif t == rdkit.Chem.rdchem.BondType.TRIPLE:
#         return '#'
#     else:
#         assert("molenc_frag.py: char_of_bond_type" == "")

# "#atoms:15 NCGC00261552-01_f00"
def read_atoms_header_line(line):
    (atoms, nb_atoms, frag_name) = [t(s) for t,s in
                                    zip((str,int,str), re.split('[: ]', line))]
    assert(atoms == "#atoms")
    return (nb_atoms, frag_name)

# "0 0,6,2,0"
def read_atom_line(line):
    (index, nb_pi, atomic_num, nb_HA, charge) = [t(s) for t,s in
                                                 zip((int,int,int,int,int),
                                                     re.split('[, ]', line))]
    return (index, nb_pi, atomic_num, nb_HA, charge)

# "#bonds:16"
def read_bonds_header_line(line):
    (bonds, nb_bonds) = [t(s) for t,s in
                         zip((str,int), re.split('[:]', line))]
    assert(bonds == "#bonds")
    return nb_bonds

# "0 - 16"
def read_bond_line(line):
    (start_i, bond_order_char, stop_i) = [t(s) for t,s in
                                          zip((int,str,int),
                                              re.split('[ ]', line))]
    return (start_i, bond_order_char, stop_i)

# "#anchors:1"
def read_anchors_header_line(line):
    (anchors, nb_anchors) = [t(s) for t,s in
                             zip((int), re.split('[:]', line))]
    assert(anchors == "#anchors")
    return nb_anchors

# "0,6,2,0 0 0,6,2,0"
def read_anchor_line(line):
    (start_t, start_i, stop_t) = [t(s) for t,s in
                                  zip((str,int,str), re.split('[ ]', line))]
    return start_i

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
