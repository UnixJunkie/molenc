#!/usr/bin/env python3

# wildcard atom scan of molecules
# for each molecule, output all variants where
# a single heavy atom at a time is switched to the SMILES wildcard atom '*'

import rdkit, sys, time
from rdkit import Chem
# from rdkit.Chem import rdChemReactions

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = " ".join(words[1:]) # everything after the SMILES string
            yield (name, smile)

if __name__ == '__main__':
    before = time.time()
    argc = len(sys.argv)
    if argc != 2:
        print("usage: %s input.smi" % sys.argv[0])
        sys.exit(1)
    input = sys.argv[1]
    count = 0
    wildcard = Chem.Atom(0)
    for name, orig_smile in RobustSmilesMolSupplier(input):
        mol = Chem.MolFromSmiles(orig_smile)
        # output original molecule first
        print("%s\t%s" % (orig_smile, name))
        num_atoms = mol.GetNumAtoms()
        # then output its variants
        for i in range(num_atoms):
            editable = Chem.EditableMol(mol)
            editable.ReplaceAtom(i, wildcard, preserveProps=True)
            edited = editable.GetMol()
            smi = Chem.MolToSmiles(edited)
            print("%s\t%s_%d" % (smi, name, i))
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)

# # original code by @Iwatobipen
# # replace any aromatic carbon to aromatic nitrogen.
# # TODO: does not compile
# def nitrogen_scan(mol_in):
#     out_mol_list = []
#     used = set()
#     rxn = rdChemReactions.ReactionFromSmarts("[c:1][H]>>[n:1]")
#     products = rxn.RunReactants([mol_in])
#     for p in products:
#         smi = Chem.MolToSmiles(Chem.RemoveHs(p))
#         if smi not in used:
#             used.add(smi)
#             out_mol_list.append(p)
#     return out_mol_list
