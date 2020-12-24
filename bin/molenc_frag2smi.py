#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# txt fragment to smi

import argparse, rdkit, re, sys, time
import molenc_common as common
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# "#anchors:1"
def read_anchors_header(line):
    (anchors, nb_anchors) = [t(s) for t,s in
                             zip((str,int), re.split('[:]', line))]
    assert(anchors == "#anchors")
    return nb_anchors

# "0,6,2,0 0 0,6,2,0"
def read_anchor(line):
    (start_t, start_i, stop_t) = [t(s) for t,s in
                                  zip((str,int,str), re.split('[ ]', line))]
    return start_i

# create a fake molecule for the corresp. fragment
def read_one_fragment(input):
    res_mol = Chem.RWMol()
    atoms_header = input.readline().strip()
    if atoms_header == '':
        raise common.End_of_file # no EOF in Python...
    nb_atoms, frag_name = common.read_atoms_header(atoms_header)
    old2new = {}
    pi_electrons = {}
    for _i in range(nb_atoms):
        line = input.readline().strip()
        (index, nb_pi, atomic_num, nb_HA, charge) = common.read_atom(line)
        # add atom
        a = Chem.Atom(atomic_num)
        a.SetFormalCharge(charge)
        j = res_mol.AddAtom(a)
        # # WARNING: we need a fix for aromatic Nitrogen ???
        # if atomic_num == 7 and nb_pi == 1 and nb_HA == 2 and charge == 0:
        #     a.SetNumExplicitHs(0)
        # we will need to convert atom indexes later
        old2new[index] = j
        # need to remember pi electrons count
        pi_electrons[j] = nb_pi
    bonds_header = input.readline().strip()
    nb_bonds = common.read_bonds_header(bonds_header)
    for i in range(nb_bonds):
        line = input.readline().strip()
        (start_i, bt, stop_i) = common.read_bond(line)
        start = old2new[start_i]
        stop = old2new[stop_i]
        # print('%d %d' % (start, stop))
        # add bond
        res_mol.AddBond(start, stop, bt)
    anchors_header = input.readline().strip()
    nb_anchors = read_anchors_header(anchors_header)
    anchors = []
    # position aromaticity flag
    for a in res_mol.GetAtoms():
        anum = a.GetAtomicNum()
        if a.IsInRing() and (anum == 6 or anum == 7):
            j = a.GetIdx()
            if pi_electrons[j] == 1:
                a.SetIsAromatic(True)
    for _i in range(nb_anchors):
        line = input.readline().strip()
        anchor = read_anchor(line)
        start = old2new[anchor]
        # dandling attachment point: dummy atom
        a = Chem.Atom('*')
        j = res_mol.AddAtom(a)
        res_mol.AddBond(start, j, Chem.rdchem.BondType.SINGLE)
        anchors.append(anchor)
    ## debug log
    # print('%s %d %d %d' % (frag_name, nb_atoms, nb_bonds, nb_anchors),
    #       file=sys.stderr)
    # smi for mol
    error = False
    try:
        Chem.SanitizeMol(res_mol)
    except rdkit.Chem.rdchem.AtomKekulizeException:
        print("AtomKekulizeException in %s" % frag_name, file=sys.stderr)
        error = True
    except rdkit.Chem.rdchem.KekulizeException:
        print("KekulizeException in %s" % frag_name, file=sys.stderr)
        error = True
    # if error:
    #     d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    #     d.drawOptions().addAtomIndices = True
    #     d.DrawMolecule(res_mol)
    #     d.FinishDrawing()
    #     png_fn = '%s.png' % frag_name
    #     with open(png_fn, 'wb') as fn:
    #         fn.write(d.GetDrawingText())
    smi = Chem.MolToSmiles(res_mol)
    return (smi, frag_name)

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
    output = open(args.output_fn, 'w')
    count = 0
    with open(input_fn) as input:
        try:
            while True:
                smi, name = read_one_fragment(input)
                count += 1
                print('%s\t%s' % (smi, name), file=output)
        except common.End_of_file:
            pass
    after = time.time()
    dt = after - before
    print("%d fragments at %.2f frag/s" % (count, count / dt), file=sys.stderr)
    output.close()
