#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# txt molecule to SMILES

import argparse, rdkit, re, sys, time
import molenc_common as common
from rdkit import Chem

# create a fake molecule for the corresp. fragment
def read_one_molecule(input):
    res_mol = Chem.RWMol()
    atoms_header = input.readline().strip()
    if atoms_header == '':
        raise common.End_of_file # no EOF in Python...
    nb_atoms, name = common.read_atoms_header(atoms_header)
    old2new = {}
    for _i in range(nb_atoms):
        line = input.readline().strip()
        (index, nb_pi, atomic_num, nb_HA, charge, stereo) = \
          common.read_atom(line)
        # add atom
        a = Chem.Atom(atomic_num)
        a.SetFormalCharge(charge)
        if stereo > 0: # set chirality
            a.SetChiralTag(common.atom_stereo_code_to_chiral_tag(stereo))
        j = res_mol.AddAtom(a)
        # we need to convert atom indexes
        old2new[index] = j
    bonds_header = input.readline().strip()
    nb_bonds = common.read_bonds_header(bonds_header)
    stereo_bonds = []
    for i in range(nb_bonds):
        line = input.readline().strip()
        (start_i, bt, stop_i, (stereo, c, d)) = common.read_bond(line)
        start = old2new[start_i]
        stop = old2new[stop_i]
        # add bond
        n = res_mol.AddBond(start, stop, bt)
        if stereo != rdkit.Chem.rdchem.BondStereo.STEREONONE:
            bi = n - 1
            # convert stereo bond stereo atoms indexes
            a = old2new[c]
            b = old2new[d]
            stereo_bonds.append((bi, stereo, a, b))
    # all atoms and bonds are here now
    # so stereo bonds info can be set
    for (bi, stereo, a, b) in stereo_bonds:
        bond = res_mol.GetBondWithIdx(bi)
        bond.SetStereo(stereo)
        bond.SetStereoAtoms(a, b)
        print('%s stereo %s on bond %d (%d, %d)' %
              (name, common.char_of_bond_stereo(stereo), bi, a, b))
    Chem.AssignStereochemistry(res_mol) # ! MANDATORY !
    try:
        Chem.SanitizeMol(res_mol)
    except rdkit.Chem.rdchem.KekulizeException:
        print("KekulizeException in %s" % name, file=sys.stderr)
    smi = Chem.MolToSmiles(res_mol)
    return (smi, name)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "txt molecule to smi")
    parser.add_argument("-i", metavar = "input.mols", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    count = 0
    with open(input_fn) as input:
        try:
            while True:
                smi, name = read_one_molecule(input)
                count += 1
                print('%s\t%s' % (smi, name), file=output)
        except common.End_of_file:
            pass
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f molecule/s" %
          (count, count / dt), file=sys.stderr)
    output.close()
