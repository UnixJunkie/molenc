#!/usr/bin/env python

# type atoms of a molecule a la atom pairs
# (nb. pi electrons if > 0, elt. symbol, nbHA neighbors)

from __future__ import print_function

import rdkit
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AtomPairs import Pairs

PeriodicTable = Chem.GetPeriodicTable()

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = words[1]
            yield (name, Chem.MolFromSmiles(smile))

def SdfMolSupplier(fn):
    for mol in Chem.SDMolSupplier(fn):
        if mol:
            name = mol.GetProp('_Name')
            yield (name, mol)

def nb_heavy_atom_neighbors(a):
    neighbors = a.GetNeighbors()
    res = 0
    for n in neighbors:
        if n.GetAtomicNum() != 1:
            res = res + 1
    return res

def type_atom(a):
    nb_pi_electrons = Pairs.Utils.NumPiElectrons(a)
    symbol = PeriodicTable.GetElementSymbol(a.GetAtomicNum())
    nbHA = nb_heavy_atom_neighbors(a)
    res = ""
    if nb_pi_electrons > 0:
        res = "%d%s%d" % (nb_pi_electrons, symbol, nbHA)
    else:
        res = "%s%d" % (symbol, nbHA)
    return res

def encode_molecule(m):
    atoms = m.GetAtoms()
    res = map(type_atom, atoms)
    return res

def print_encoded_atoms(atoms):
    for i, a in enumerate(atoms):
        if i > 0:
            print(' %s' % a, end='')
        else:
            print(a, end='')

def main():
    if len(sys.argv) != 2:
        print("usage: %s input.{smi|sdf}" % sys.argv[0])
        sys.exit(1)
    input = sys.argv[1]
    mol_supplier = None
    if input.endswith(".smi"):
        mol_supplier = RobustSmilesMolSupplier
    if input.endswith(".sdf"):
        mol_supplier = SdfMolSupplier
    for name, mol in mol_supplier(input):
        if mol is None:
            continue
        print_encoded_atoms(encode_molecule(mol))
        print('\t%s' % name)

if __name__ == '__main__':
    main()
