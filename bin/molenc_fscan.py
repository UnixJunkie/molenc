#!/usr/bin/env python3
#
# Fluorine scan of a molecule
# create all analogs of input molecule where one heavy atom
# at a time has all of its hydrogens replaced by F
# optionally, do this only for heteroatoms

import rdkit, sys
from rdkit import Chem

def RobustSmilesMolSupplier(input_fn):
    with open(input_fn) as f:
        for line in f:
            strip = line.strip()
            toks = strip.split()
            smi = toks[0]
            toks.reverse()
            name = toks[0]
            yield (smi, name)

fluor = Chem.Atom(9)

# FBR: proper CLI handling
if __name__ == '__main__':
    argc = len(sys.argv)
    if argc != 2:
        print("usage: %s input.smi [--hetero]" % sys.argv[0])
        sys.exit(1)
    input_fn = sys.argv[1]
    for smi, name in RobustSmilesMolSupplier(input_fn):
        # output original molecule first
        print("%s\t%s" % (smi, name))
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        # then output its variants
        count = 1
        for a in mol.GetAtoms():
            if a.GetAtomicNum() > 1:
                # heavy atom
                if a.GetTotalNumHs(includeNeighbors=True) >= 1:
                    # hydrogens attached
                    editable = Chem.EditableMol(mol)
                    for neighb in a.GetNeighbors():
                        if neighb.GetAtomicNum() == 1:
                            # Fluorine instead
                            a_j = neighb.GetIdx()
                            editable.ReplaceAtom(a_j, fluor)
                    edited = editable.GetMol()
                    smi = Chem.MolToSmiles(edited)
                    print("%s\t%s_%d" % (smi, name, count))
                    count += 1
