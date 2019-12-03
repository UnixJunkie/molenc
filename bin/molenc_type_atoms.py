#!/usr/bin/env python3

# type atoms of a molecule a la atom pairs
# (nb. pi electrons if > 0, elt. symbol, nbHA neighbors)
# formal charges are ignored, as was the case in the seminal implementation
# of atom pairs, not sure this is very smart though

import molenc_common, os, rdkit, sys, time
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

def SdfMolSupplier(fn):
    for mol in Chem.SDMolSupplier(fn):
        if mol:
            name = mol.GetProp('_Name')
            yield (name, mol)

def nb_heavy_atom_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() != 1:
            res += 1
    return res

PeriodicTable = Chem.GetPeriodicTable()

def string_of_charge(charge):
    if charge == 0: return ""
    elif charge == -1: return "-"
    elif charge == 1: return "+"
    else: return ("%+d" % charge)

def type_atom(a):
    res = None
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    symbol = PeriodicTable.GetElementSymbol(a.GetAtomicNum())
    nbHA = nb_heavy_atom_neighbors(a)
    formal_charge = string_of_charge(a.GetFormalCharge())
    if nb_pi_electrons > 0:
        res = "%d%s%d%s" % (nb_pi_electrons, symbol, nbHA, formal_charge)
    else:
        res = "%s%d%s" % (symbol, nbHA, formal_charge)
    return res

# stereo information encoding, at the atom level and using integers
class AtomStereoCode(Enum):
    NONE = 0
    ANY_CENTER = 1
    R_CENTER = 2
    S_CENTER = 3
    ANY_BOND_END = 4
    Z_BOND_END = 5
    E_BOND_END = 6
    CIS_BOND_END = 7
    TRANS_BOND_END = 8

to_atom_stereo = { '?': AtomStereoCode.ANY_CENTER,
                   'R': AtomStereoCode.R_CENTER,
                   'S': AtomStereoCode.S_CENTER,
                   rdkit.Chem.rdchem.BondStereo.STEREONONE: AtomStereoCode.NONE,
                   rdkit.Chem.rdchem.BondStereo.STEREOANY: AtomStereoCode.ANY_BOND_END,
                   rdkit.Chem.rdchem.BondStereo.STEREOZ: AtomStereoCode.Z_BOND_END,
                   rdkit.Chem.rdchem.BondStereo.STEREOE: AtomStereoCode.E_BOND_END,
                   rdkit.Chem.rdchem.BondStereo.STEREOCIS: AtomStereoCode.CIS_BOND_END,
                   rdkit.Chem.rdchem.BondStereo.STEREOTRANS: AtomStereoCode.TRANS_BOND_END }

def get_stereo_codes(m):
    # by default, each atom has no stereo
    res = [AtomStereoCode.NONE for i in range(m.GetNumAtoms())]
    # unless detected otherwise for stereo bonds
    for b in m.GetBonds():
        bstereo = b.GetStereo()
        if bstereo != rdkit.Chem.rdchem.BondStereo.STEREONONE:
            i = b.GetBeginAtomIdx()
            j = b.GetEndAtomIdx()
            k = to_atom_stereo[bstereo]
            res[i] = k
            res[j] = k
    # or for chiral centers
    for i, k in Chem.FindMolChiralCenters(m):
        res[i] = to_atom_stereo[k]
    return res

# # tests
# cis = Chem.MolFromSmiles('C/C=C\C')
# trans = Chem.MolFromSmiles('C/C=C/C')
# l_ala = Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
# d_ala = Chem.MolFromSmiles('N[C@H](C)C(=O)O')
# print(get_stereo_codes(cis))
# print(get_stereo_codes(trans))
# print(get_stereo_codes(l_ala))
# print(get_stereo_codes(d_ala))

def encode_molecule(m):
    return map(type_atom, m.GetAtoms())

def print_encoded_atoms(atoms):
    for i, a in enumerate(atoms):
        print("%d %s" % (i, a))

if __name__ == '__main__':
    before = time.time()
    argc = len(sys.argv)
    if argc != 2:
        print("usage: %s input.{smi|sdf}" % sys.argv[0])
        sys.exit(1)
    input = sys.argv[1]
    mol_supplier = None
    if input.endswith(".smi"):
        mol_supplier = RobustSmilesMolSupplier
    if input.endswith(".sdf"):
        mol_supplier = SdfMolSupplier
    count = 0
    for name, mol in mol_supplier(input):
        print("#atoms:%d %s" % (mol.GetNumAtoms(), name))
        print_encoded_atoms(encode_molecule(mol))
        molenc_common.print_bonds(mol)
        molenc_common.print_distance_matrix(mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
