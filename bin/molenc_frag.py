#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# atom typing and molecule fragmentation hints

import argparse
import molenc_common as common
import rdkit
import sys
import time

from common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from common import StereoCodes

to_stereo_code = \
  { # atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    '?': StereoCodes.ANY_CENTER,
    # SMILES @@ means clockwise / R / Chem.ChiralType.CHI_TETRAHEDRAL_CW
    'R': StereoCodes.R_CENTER,
    # SMILES @ means anti clockwise / S / Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    'S': StereoCodes.S_CENTER,
    rdkit.Chem.rdchem.BondStereo.STEREONONE: StereoCodes.NONE,
    rdkit.Chem.rdchem.BondStereo.STEREOANY: StereoCodes.ANY_BOND,
    rdkit.Chem.rdchem.BondStereo.STEREOZ: StereoCodes.Z_BOND,
    rdkit.Chem.rdchem.BondStereo.STEREOE: StereoCodes.E_BOND,
    rdkit.Chem.rdchem.BondStereo.STEREOCIS: StereoCodes.CIS_BOND,
    rdkit.Chem.rdchem.BondStereo.STEREOTRANS: StereoCodes.TRANS_BOND }

def get_atom_stereo_codes(m):
    # by default, each atom has no stereo
    res = [StereoCodes.NONE for i in range(m.GetNumAtoms())]
    # # unless detected otherwise for stereo bonds
    # for b in m.GetBonds():
    #     bstereo = b.GetStereo()
    #     if bstereo != rdkit.Chem.rdchem.BondStereo.STEREONONE:
    #         i = b.GetBeginAtomIdx()
    #         j = b.GetEndAtomIdx()
    #         k = to_stereo_code[bstereo]
    #         res[i] = k
    #         res[j] = k
    # or for chiral centers
    for i, k in Chem.FindMolChiralCenters(m):
        res[i] = to_stereo_code[k]
    return res

# # stereo code tests
# cis = Chem.MolFromSmiles('C/C=C\C')
# trans = Chem.MolFromSmiles('C/C=C/C')
# l_ala = Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
# d_ala = Chem.MolFromSmiles('N[C@H](C)C(=O)O')
# print(get_stereo_codes(cis))
# print(get_stereo_codes(trans))
# print(get_stereo_codes(l_ala))
# print(get_stereo_codes(d_ala))

def print_typed_atoms(out, mol):
    stereo = get_atom_stereo_codes(mol)
    for a in mol.GetAtoms():
        i = a.GetIdx()
        t = common.type_atom(a)
        s = stereo[i]
        print("%d %s,%d" % (i, t, s), file=out)

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

def string_of_bond_stereo(bond):
    st = bond.GetStereo()
    c = common.char_of_bond_stereo(st)
    if c == 'N':
        return c
    else:
        (a, b) = bond.GetStereoAtoms()
        str = "%c:%d:%d" % (c, a, b)
        return str

# print all bonds with their type (and optional stereo info)
def print_bonds(out, mol):
    print("#bonds:%d" % mol.GetNumBonds(), file=out)
    bonds = mol.GetBonds()
    for bond in bonds:
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        t = char_of_bond_type(bond)
        stereo = string_of_bond_stereo(bond)
        print("%d %c %d %s" % (a, t, b, stereo), file=out)

# print which bonds are cuttable and the suggested number of cuts
def print_cuttable_bonds(out, mol):
    cuttable_bonds = common.find_cuttable_bonds(mol)
    total_weight = Descriptors.MolWt(mol)
    # 150 Da: D. Rognan's suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cuttable_bonds), nb_frags - 1)
    print("#cut_bonds:%d:%d" % (len(cuttable_bonds), max_cuts), file=out)
    for bond in cuttable_bonds:
        i = bond.GetIdx()
        print("%d" % i, file=out)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "compute molecule fragmentation hints")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.txt", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--draw", dest = "draw_mol", action ='store_true',
                        default = False,
                        help = "output PNG for each molecule w/ atom indexes")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    draw_mol = args.draw_mol
    output = open(args.output_fn, 'w')
    count = 0
    # fragmenting ---------------------------------------------------------
    mol_supplier = RobustSmilesMolSupplier(input_fn)
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
