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
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# "#anchors:1"
def read_anchors_header(line):
    (anchors, nb_anchors) = [t(s) for t,s in
                             zip((str,int), re.split('[:]', line))]
    assert(anchors == "#anchors")
    return nb_anchors

# "0,6,2,0,0 0 0,6,2,0,0"
def read_anchor(line):
    (start_t, start_i, stop_t) = [t(s) for t,s in
                                  zip((str,int,str), re.split('[ ]', line))]
    return start_i

def current_valence(a):
    res = 0.0
    for b in a.GetBonds():
        bt = b.GetBondType()
        if bt == rdkit.Chem.rdchem.BondType.SINGLE:
            res += 1.0
        elif bt == rdkit.Chem.rdchem.BondType.AROMATIC:
            res += 1.5
        elif bt == rdkit.Chem.rdchem.BondType.DOUBLE:
            res += 2.0
        elif bt == rdkit.Chem.rdchem.BondType.TRIPLE:
            res += 3.0
    return res

def count_aromatic_bonds(a):
    count = 0
    for b in a.GetBonds():
        if b.GetBondType() == rdkit.Chem.rdchem.BondType.AROMATIC:
            count += 1
    return count

# try to correct a molecule that failed kekulization
def restore_from_kekulization_error(mol):
    print('---------', file=sys.stderr)
    pat1 = Chem.MolFromSmarts('c(=O)')
    rep1 = Chem.MolFromSmarts('c(O)')
    res_mol1 = AllChem.ReplaceSubstructs(mol, pat1, rep1)
    res_mol = res_mol1[0]
    for a in res_mol.GetAtoms():
        anum = a.GetAtomicNum()
        val = current_valence(a)
        if anum == 7 and (not a.IsInRing()):
            print('N not in ring with %d aromatic bonds; val=%.1f' %
                  (count_aromatic_bonds(a), val), file = sys.stderr)
        if a.IsInRing() and anum == 7:
            print('N in ring with %d aromatic bonds; val=%.1f' %
                  (count_aromatic_bonds(a), val), file = sys.stderr)
            if val == 4.0:
                # set the correct +1 partial charge
                # for Nitrogen with valence 4 in rings
                a.SetFormalCharge(1)
                # forbid addition of implicit hydrogens later
                a.SetNoImplicit(True)
                a.SetNumExplicitHs(0)
                print('N charge correct', file = sys.stderr)
            # if count_aromatic_bonds(a) == 2 and val == 3.0:
            #     # aromatic Nitrogen _must_ be [nH] instead of [n]
            #     a.SetNumExplicitHs(1)
            if count_aromatic_bonds(a) >= 1:
                # aromatic Nitrogen _must_ be [nH] instead of [n]
                a.SetNumExplicitHs(1)
    # pat2 = Chem.MolFromSmarts('[n+0]')
    # rep2 = Chem.MolFromSmarts('[nH]')
    # res_mol2 = AllChem.ReplaceSubstructs(res_mol, pat2, rep2)
    # res_mol = res_mol2[0]
    return res_mol

# create a fake molecule for the corresp. fragment
def read_one_fragment(input):
    res_mol = Chem.RWMol()
    atoms_header = input.readline().strip()
    if atoms_header == '':
        raise common.End_of_file # no EOF in Python...
    nb_atoms, frag_name = common.read_atoms_header(atoms_header)
    old2new = {}
    # read atoms
    for _i in range(nb_atoms):
        line = input.readline().strip()
        (index, nb_pi, atomic_num, nb_HA, charge, stereo) = common.read_atom(line)
        # add atom
        a = Chem.Atom(atomic_num)
        a.SetFormalCharge(charge)
        if stereo > 0: # set chirality
            a.SetChiralTag(common.atom_stereo_code_to_chiral_tag(stereo))
        j = res_mol.AddAtom(a)
        # we will need to convert atom indexes later
        old2new[index] = j
    bonds_header = input.readline().strip()
    nb_bonds = common.read_bonds_header(bonds_header)
    stereo_bonds = []
    # PREVENT CUTTING BONDS BETWEEN A STEREO BOND AND ONE
    # OF ITS STEREO ATOMS DURING FRAGMENTATION !!!
    # read bonds
    for i in range(nb_bonds):
        line = input.readline().strip()
        (start_i, bt, stop_i, (stereo, c, d)) = common.read_bond(line)
        # convert atom indexes
        start = old2new[start_i]
        stop = old2new[stop_i]
        # print('%d %d' % (start, stop))
        # add bond
        n = res_mol.AddBond(start, stop, bt)
        if stereo != rdkit.Chem.rdchem.BondStereo.STEREONONE:
            bi = n - 1
            # if an exception is thrown here, it means
            # this stereobond stereo atom is out of the fragments atoms
            # (attachment points dummy atoms are forbidden)
            #
            # convert stereo bond stereo atoms indexes
            a = old2new[c]
            b = old2new[d]
            stereo_bonds.append((bi, stereo, a, b))
    # all fragments atoms and internal bonds are here now
    # so stereo bonds info can be set
    for (bi, stereo, a, b) in stereo_bonds:
        bond = res_mol.GetBondWithIdx(bi)
        bond.SetStereo(stereo)
        bond.SetStereoAtoms(a, b)
        #print('%s stereo on bond %d (%d, %d)' % (frag_name, bi, a, b))
    anchors_header = input.readline().strip()
    nb_anchors = read_anchors_header(anchors_header)
    anchors = []
    for _i in range(nb_anchors):
        line = input.readline().strip()
        anchor = read_anchor(line)
        start = old2new[anchor]
        # dandling attachment point: dummy atom
        a = Chem.Atom('*')
        j = res_mol.AddAtom(a)
        # only single, non stereo, bonds out of rings have been cut
        res_mol.AddBond(start, j, Chem.rdchem.BondType.SINGLE)
        anchors.append(anchor)
    # smi for mol
    try:
        Chem.SanitizeMol(res_mol)
        res_smi = Chem.MolToSmiles(res_mol)
        return (False, res_smi, frag_name)
    except rdkit.Chem.rdchem.KekulizeException:
        print("KekulizeException in %s" % frag_name, file=sys.stderr)
        return (True, "", frag_name)
    # if kekul_error:
    #     res_mol = restore_from_kekulization_error(res_mol)
    #     try:
    #         Chem.SanitizeMol(res_mol)
    #         # print('restored one')
    #     except rdkit.Chem.rdchem.KekulizeException:
    #         print('kekulization rescue failed\n%s\t%s' %
    #               (Chem.MolToSmiles(res_mol), frag_name), file=sys.stderr)
    #         kekul_error2 = True
    #     ## draw mol with atom indexes for inspection
    #     #   d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    #     #   d.drawOptions().addAtomIndices = True
    #     #   d.DrawMolecule(res_mol)
    #     #   d.FinishDrawing()
    #     #   png_fn = '%s.png' % frag_name
    #     #   with open(png_fn, 'wb') as fn:
    #     #     fn.write(d.GetDrawingText())

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = "convert txt fragments to SMILES")
    parser.add_argument("-i", metavar = "input.frags", dest = "input_fn",
                        help = "fragments input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "SMILES output file")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of how to use
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output = open(args.output_fn, 'w')
    err_out = open(args.output_fn + '.err.names', 'w')
    count = 0
    kekul_err_count = 0
    with open(input_fn) as input:
        try:
            while True:
                kekul_err, smi, name = read_one_fragment(input)
                count += 1
                if kekul_err:
                    print('%s' % name, file=err_out)
                    kekul_err_count += 1
                else:
                    print('%s\t%s' % (smi, name), file=output)
        except common.End_of_file:
            pass
    after = time.time()
    dt = after - before
    print("%d (%d kekul. errors) fragments at %.2f frag/s" %
          (count, kekul_err_count, count / dt), file=sys.stderr)
    output.close()
    err_out.close()
