#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.
#
# Compute the Bemis-Murcho generic scaffold (framework)
# of each input molecule.
#
# Bemis, G. W., & Murcko, M. A. (1996).
# "The properties of known drugs. 1. Molecular frameworks."
# Journal of medicinal chemistry, 39(15), 2887-2893.

import argparse, rdkit, sys
from rdkit import Chem

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smi = words[0]
            name = words[1]
            mol = Chem.MolFromSmiles(smi)
            yield (smi, name, mol)

def find_terminal_atoms(mol):
    res = []
    for a in mol.GetAtoms():
        if len(a.GetBonds()) == 1:
            res.append(a)
    return res

def BemisMurckoFramework(mol):
    # keep only Heavy Atoms (HA)
    only_HA = rdkit.Chem.rdmolops.RemoveHs(mol)
    # switch all HA to Carbon
    rw_mol = Chem.RWMol(only_HA)
    for i in range(rw_mol.GetNumAtoms()):
        rw_mol.ReplaceAtom(i, Chem.Atom(6))
    # switch all non single bonds to single
    non_single_bonds = []
    for b in rw_mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE:
            non_single_bonds.append(b)
    for b in non_single_bonds:
        j = b.GetBeginAtomIdx()
        k = b.GetEndAtomIdx()
        rw_mol.RemoveBond(j, k)
        rw_mol.AddBond(j, k, Chem.BondType.SINGLE)
    # as long as there are terminal atoms, remove them
    terminal_atoms = find_terminal_atoms(rw_mol)
    while terminal_atoms != []:
        for a in terminal_atoms:
            for b in a.GetBonds():
                rw_mol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            rw_mol.RemoveAtom(a.GetIdx())
        terminal_atoms = find_terminal_atoms(rw_mol)
    return rw_mol.GetMol()

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Append Bemis-Murcko scaffold to each input molecule")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file")
    parser.add_argument("-o", metavar = "output_smi", dest = "output_smi",
                        help = "output SMILES file")
    parser.add_argument('--new-line', dest='new_line',
                        action='store_true', default=False,
                        help = "insert a newline before the scaffold")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_smi = args.output_smi
    new_line = args.new_line
    out_count = 0
    error_count = 0
    with open(output_smi, 'w') as out_file:
        for smi, name, mol in RobustSmilesMolSupplier(input_smi):
            if mol is None:
                error_count += 1
            else:
                scaff = BemisMurckoFramework(mol)
                scaff_smi = Chem.MolToSmiles(scaff)
                if new_line:
                    print("%s\t%s\n%s" % (smi, name, scaff_smi), file=out_file)
                else:
                    print("%s\t%s\t%s" % (smi, name, scaff_smi), file=out_file)
                out_count += 1
    total_count = out_count + error_count
    print("encoded: %d errors: %d total: %d" %
          (out_count, error_count, total_count),
          file=sys.stderr)

if __name__ == '__main__':
    main()
