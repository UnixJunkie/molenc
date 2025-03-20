#!/usr/bin/env python3
#
# Fluorine scan of a molecule
# create all analogs of input molecule where one heavy atom
# at a time has all of its hydrogens replaced by F
# optionally, do this only for heteroatoms

import argparse, rdkit, sys
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

if __name__ == '__main__':
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "compute atom types and distances")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "molecules output file")
    parser.add_argument('--hetero', dest='only_heteroatoms', action='store_true',
                        help = "only scan heteroatoms")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    only_hetero = args.only_heteroatoms
    with open(output_fn, 'w') as out:
        for smi, name in RobustSmilesMolSupplier(input_fn):
            # output original molecule first
            print("%s\t%s" % (smi, name), file=out)
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            # then output its fluorinated analogs
            count = 1
            for a in mol.GetAtoms():
                anum = a.GetAtomicNum()
                if anum > 1 and ((not only_hetero) or anum != 6):
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
                        print("%s\t%s_%d" % (smi, name, count), file=out)
                        count += 1
