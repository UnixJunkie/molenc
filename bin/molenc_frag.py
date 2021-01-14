#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# atom typing and molecule fragmentation hints

import argparse, molenc_common, os, random, rdkit, sys, time

from enum import IntEnum
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Draw import rdMolDraw2D
from molenc_common import StereoCodes

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = " ".join(words[1:]) # everything after the SMILES string
            mol = Chem.MolFromSmiles(smile)
            mol.SetProp("name", name)
            yield (name, mol)

def nb_heavy_atom_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() != 1:
            res += 1
    return res

to_stereo_code = \
  { '?': StereoCodes.ANY_CENTER, # atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    'R': StereoCodes.R_CENTER, # SMILES @@ means clockwise / R / Chem.ChiralType.CHI_TETRAHEDRAL_CW
    'S': StereoCodes.S_CENTER, # SMILES @ means anti clockwise / S / Chem.ChiralType.CHI_TETRAHEDRAL_CCW
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

def type_atom(a):
    # stereo chemistry is ignored for the moment
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    atom_num = a.GetAtomicNum()
    nbHA = nb_heavy_atom_neighbors(a)
    formal_charge = a.GetFormalCharge()
    # make this easy to parse / unambiguous
    res = "%d,%d,%d,%d" % (nb_pi_electrons, atom_num, nbHA, formal_charge)
    return res

def log_protected_bond(name, b):
    print('mol %s: protected bond %d' % (name, b.GetIdx()))

# no stereo, single bonds not in rings
def find_cuttable_bonds(mol):
    # protect bonds between stereo bond atoms and their stereo atoms
    name = mol.GetProp("name")
    for b in mol.GetBonds():
        if b.GetStereo() != rdkit.Chem.rdchem.BondStereo.STEREONONE:
            (i, j) = (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            (k, l) = b.GetStereoAtoms()
            b0 = mol.GetBondBetweenAtoms(i, k)
            b1 = mol.GetBondBetweenAtoms(i, l)
            b2 = mol.GetBondBetweenAtoms(j, k)
            b3 = mol.GetBondBetweenAtoms(j, l)
            if b0 != None:
                b0.SetBoolProp("protected", True)
                log_protected_bond(name, b0)
            if b1 != None:
                b1.SetBoolProp("protected", True)
                log_protected_bond(name, b1)
            if b2 != None:
                b2.SetBoolProp("protected", True)
                log_protected_bond(name, b2)
            if b3 != None:
                b3.SetBoolProp("protected", True)
                log_protected_bond(name, b3)
    res = []
    for b in mol.GetBonds():
        if ((b.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE) and
            (not b.IsInRing()) and
            (b.GetStereo() == rdkit.Chem.rdchem.BondStereo.STEREONONE) and
            (b.HasProp("protected") == 0)): # HasProp returns an int... :(
            res.append(b)
    return res

def print_typed_atoms(out, mol):
    stereo = get_atom_stereo_codes(mol)
    for a in mol.GetAtoms():
        i = a.GetIdx()
        t = type_atom(a)
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
    c = molenc_common.char_of_bond_stereo(st)
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
    cuttable_bonds = find_cuttable_bonds(mol)
    total_weight = Descriptors.MolWt(mol)
    # 150 Da: D. Rognan's suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cuttable_bonds), nb_frags)
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
