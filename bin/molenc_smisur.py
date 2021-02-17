#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# The Smiling Surgeon: a doctor operating directly at the SMILES level

import argparse
import ast
import molenc_common as common
import random
import rdkit
import sys
import time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import RWMol
from rdkit.Chem.AtomPairs import Pairs

def get_name(mol):
    return mol.GetProp("name")

def set_name(mol, name):
    mol.SetProp("name", name)

def index_for_atom_type(atom_types_dict, atom_type):
    try:
        return atom_types_dict[atom_type]
    except KeyError:
        # want indexes to start at 1; so the isotope number is
        # always explicit in the fragments SMILES output
        v = len(atom_types_dict) + 1
        atom_types_dict[atom_type] = v
        return v

def dict_reverse_binding(dico):
    res = {}
    for k, v in dico.items():
        res[v] = k
    return res

def fragment_on_bonds_and_label(mol, bonds):
    labels = []
    atom_type_to_index = {}
    for bi in bonds:
        b = mol.GetBondWithIdx(bi)
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        # get or create dictionary keys for those atom types
        ai = mol.GetAtomWithIdx(i)
        aj = mol.GetAtomWithIdx(j)
        at_i = common.type_atom(ai)
        at_j = common.type_atom(aj)
        vi = index_for_atom_type(atom_type_to_index, at_i)
        vj = index_for_atom_type(atom_type_to_index, at_j)
        labels.append((vi, vj))
    fragmented = Chem.FragmentOnBonds(mol, bonds, dummyLabels=labels)
    smi = Chem.MolToSmiles(fragmented)
    name = get_name(mol)
    index_to_atom_type = dict_reverse_binding(atom_type_to_index)
    return (smi, name, index_to_atom_type)

# SMILES fragmentation
def cut_some_bonds(mol, seed):
    cuttable_bonds = common.find_cuttable_bonds(mol)
    cut_bonds_indexes = [b.GetIdx() for b in cuttable_bonds]
    total_weight = Descriptors.MolWt(mol)
    # 150 Da: D. Rognan's suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cut_bonds_indexes), nb_frags - 1)
    # print("mol %s; cut %d bonds" % (mol.GetProp("name"), max_cuts),
    #       file=sys.stderr)
    random.shuffle(cut_bonds_indexes)
    to_cut = cut_bonds_indexes[0:max_cuts]
    if len(to_cut) == 0:
        # molecule too small: not fragmented
        # still, we output it so that input and output SMILES files can be
        # visualized side-by-side
        smi = Chem.MolToSmiles(mol)
        name = get_name(mol)
        dico = {}
        return (smi, name, dico)
    else:
        return fragment_on_bonds_and_label(mol, to_cut)

def FragmentsSupplier(filename):
    with open(filename) as f:
        for line in f:
            mixture, name = line.strip().split("\t") # enforce TAB-separated
            fragments_smiles = mixture.split(".")
            parent_mol_name, dict_str = name.split(";")
            dico = ast.literal_eval(dict_str)
            for i, smi in enumerate(fragments_smiles):
                frag_name = '%s_f%d' % (parent_mol_name, i)
                if len(dico) > 0: # molecule _was_ fragmented (not too small)
                    # print(smi, frag_name, dico, file=sys.stderr) # debug
                    yield (smi, frag_name, dico)

def read_all_fragments(fn):
    return [(smi, name, dico) for (smi, name, dico) in FragmentsSupplier(fn)]

def random_choose_one(all_frags):
    n = len(all_frags)
    i = random.randint(0, n - 1)
    return all_frags[i]

def count_uniq_fragment(all_frags):
    ss = set() # string set
    for (smi, _name, _dico) in all_frags:
        ss.add(smi)
    return len(ss)

# the only one attached to a dummy atom / attachment point
def get_src_atom_idx(a):
    neighbs = a.GetNeighbors()
    assert(len(neighbs) == 1)
    src_a = neighbs[0]
    return src_a.GetIdx()

# set "name" prop. to frag_mol
# record (frag_mol, src_atom_idx, src_atom_typ, dst_atom_typ)
def index_fragments(frags):
    res = {}
    for smi, frag_name, dico in frags:
        frag_mol = Chem.MolFromSmiles(smi)
        set_name(frag_mol, frag_name)
        # process each attachment point
        for a in frag_mol.GetAtoms():
            if a.GetAtomicNum() == 0: # '*' wildcard atom
                isotope = a.GetIsotope()
                # print(isotope, dico) # debug
                dst_typ = dico[isotope]
                dst_idx = a.GetIdx()
                src_idx = get_src_atom_idx(a)
                src_a = frag_mol.GetAtomWithIdx(src_idx)
                src_typ = common.type_atom(src_a)
                # record the fragment under key: (dst_typ, src_typ)
                # i.e. ready to use by requiring fragment
                key = (dst_typ, src_typ)
                # print('insert key: %s' % str(key)) # debug
                value = (frag_mol, dst_idx)
                a.SetProp("dst_typ", dst_typ)
                try:
                    previous_frags = res[key]
                    previous_frags.append(value)
                except KeyError:
                    res[key] = [value]
    return res

# extract fragments from values of the dictionary
def extract_fragments(dico):
    res = []
    for _k, v in dico.items():
        for (frag_mol, _dst_idx) in v:
            res.append(frag_mol)
    return res

# return a new molecule, where m1 and m2 are now attached
# via a single bond; after this bond is introduced, the former
# corresponding attachment points/atoms are removed
def bind_molecules(m1, m2):
    # print('m1: %s' % Chem.MolToSmiles(m1)) #debug
    # print('m2: %s' % Chem.MolToSmiles(m2)) #debug
    n1 = m1.GetNumAtoms()
    n2 = m2.GetNumAtoms()
    m = n1 + n2
    rw_mol = Chem.RWMol(Chem.CombineMols(m1, m2))
    assert(rw_mol.GetNumAtoms() == m)
    name1 = get_name(m1)
    name2 = get_name(m2)
    new_name = '%s,%s' % (name1, name2)
    set_name(rw_mol, new_name)
    # find first compatible attach points between m1 and m2 in rw_mol
    for i in range(n1):
        ai = rw_mol.GetAtomWithIdx(i)
        if ai.GetAtomicNum() == 0: # attachment point
            dst_idx = ai.GetIdx()
            dst_typ = ai.GetProp("dst_typ")
            src_idx = get_src_atom_idx(ai)
            src_a = rw_mol.GetAtomWithIdx(src_idx)
            src_typ = common.type_atom(src_a)
            j = n1
            while j < m:
                aj = rw_mol.GetAtomWithIdx(j)
                if aj.GetAtomicNum() == 0: # attachment point
                    dst_idx2 = aj.GetIdx()
                    dst_typ2 = aj.GetProp("dst_typ")
                    src_idx2 = get_src_atom_idx(aj)
                    src_a2 = rw_mol.GetAtomWithIdx(src_idx2)
                    src_typ2 = common.type_atom(src_a2)
                    if (dst_typ == src_typ2 and
                        dst_typ2 == src_typ):
                        # attach. points are compatible
                        rw_mol.AddBond(
                            src_idx, src_idx2, Chem.rdchem.BondType.SINGLE)
                        # remove former attachment points
                        rw_mol.RemoveAtom(dst_idx2) # to not shift lower atom indexes
                        rw_mol.RemoveAtom(dst_idx)
                        return rw_mol
                j += 1
    # crash if none was found
    print("bind_molecules: could not connect fragment %s w/ %s" %
          (name1, name2), file=sys.stderr)
    assert(False)

# first attach. point/atom index, or -1 if no more
def find_first_attach_index(mol):
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0: # '*' wildcard atom
            return a.GetIdx()
    return -1

# attach matching fragments until no attachment points are left
# WARNING: this is a recursive function
def grow_fragment(frag_seed_mol, frags_index):
    dst_idx = find_first_attach_index(frag_seed_mol)
    if dst_idx == -1:
        try:
            Chem.SanitizeMol(frag_seed_mol)
            # the constituting fragments might have some stereo info
            # that we want to preserve up to the final molecule
            Chem.AssignStereochemistry(frag_seed_mol) # ! MANDATORY _AFTER_ SanitizeMol !
            # print('after sanitize then stereo: %s' % Chem.MolToSmiles(res_mol), file=sys.stderr)
            return frag_seed_mol.GetMol()
        except rdkit.Chem.rdchem.KekulizeException:
            print("KekulizeException in %s" % get_name(frag_seed_mol), file=sys.stderr)
            return frag_seed_mol.GetMol()
    else:
        dst_a = frag_seed_mol.GetAtomWithIdx(dst_idx)
        dst_typ = dst_a.GetProp("dst_typ")
        src_idx = get_src_atom_idx(dst_a)
        src_a = frag_seed_mol.GetAtomWithIdx(src_idx)
        src_typ = common.type_atom(src_a)
        # draw compatible fragment
        key = (src_typ, dst_typ) # current to wanted direction
        # print('want key: %s' % str(key)) # debug
        possible_compat_frags = frags_index[key]
        compat_frag = random_choose_one(possible_compat_frags)
        (frag_mol2, dst_idx2) = compat_frag
        dst_a2 = frag_mol2.GetAtomWithIdx(dst_idx2)
        dst_typ2 = dst_a2.GetProp("dst_typ")
        src_idx2 = get_src_atom_idx(dst_a2)
        src_a2 = frag_mol2.GetAtomWithIdx(src_idx2)
        src_typ2 = common.type_atom(src_a2)
        # check fragments compatibility
        assert(src_typ == dst_typ2)
        assert(dst_typ == src_typ2)
        # connect them
        new_mol = bind_molecules(frag_seed_mol, frag_mol2)
        # rec. call
        return grow_fragment(new_mol, frags_index)

def write_out(gen_mol, count, gen_smi, output):
    frag_names = get_name(gen_mol)
    name_prfx = 'genmol%d' % count
    print("%s\t%s:%s" % (gen_smi, name_prfx, frag_names), file=output)

# Oprea's lead-like filter
# Hann, M. M., & Oprea, T. I. (2004).
# Pursuing the leadlikeness concept in pharmaceutical research.
# Current opinion in chemical biology, 8(3), 255-263.
def lead_like_filter(mol):
    # MolW <= 460
    if Descriptors.MolWt(mol) > 460:
        return False
    # -4.0 <= LogP <= 4.2
    LogP = Descriptors.MolLogP(mol)
    if LogP < -4.0 or LogP > 4.2:
        return False
    # # LogSw >= -5 # ignored
    # rotB <= 10
    if Descriptors.NumRotatableBonds(mol) > 10:
        return False
    # nRings <= 4 (number of SSSR rings, _not_ aromatic rings)
    if Chem.GetSSSR(mol) > 4:
        return False
    # HBD <= 5
    if Descriptors.NumHDonors(mol) > 5:
        return False
    # HBA <= 9
    if Descriptors.NumHAcceptors(mol) > 9:
        return False
    return True # lead-like then!

def new_enough(filter_diverse, gen_smi, seen_smiles):
    if not filter_diverse:
        return True
    else:
        res = not (gen_smi in seen_smiles)
        seen_smiles.add(gen_smi)
        return res

def lead_like_enough(ll_filter, mol):
    return ((not ll_filter) or lead_like_filter(mol))

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules, or assemble molecular fragments")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--seed", dest = "seed", default = -1,
                        type = int, help = "RNG seed")
    parser.add_argument("-n", dest = "nb_passes", default = 1,
                        type = int, help = "number of fragmentation passes")
    parser.add_argument("--assemble", dest = "nmols", default = -1,
                        type = int, help = "number of molecules to generate")
    parser.add_argument("--diverse", dest = "diverse", default = False,
                        action = "store_true",
                        help = "enforce uniqueness of generated SMILES")
    parser.add_argument("--lead-like", dest = "ll_filter", default = False,
                        action = "store_true",
                        help = "only generate lead-like molecules")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    nb_passes = args.nb_passes
    nmols = args.nmols
    assemble = nmols > 0
    rng_seed = args.seed
    diverse = args.diverse
    ll_filter = args.ll_filter
    seen_smiles = set()
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    output = open(args.output_fn, 'w')
    count = 0
    if assemble: # assembling fragments ---------------------------------------
        smi_fragments = read_all_fragments(input_fn)
        nb_uniq = count_uniq_fragment(smi_fragments)
        print('read %d fragments (uniq: %d)' % (len(smi_fragments), nb_uniq))
        index = index_fragments(smi_fragments)
        fragments = extract_fragments(index)
        print('%d fragment keys' % len(index))
        # # inspect the index (to debug)
        # for k, v in index.items():
        #     print("k:%s -> %d frags" % (k, len(v)))
        while count < nmols:
            seed_frag = random_choose_one(fragments)
            # print('seed_frag: %s' % get_name(seed_frag)) # debug
            gen_mol = grow_fragment(seed_frag, index)
            gen_smi = Chem.MolToSmiles(gen_mol)
            is_new = new_enough(diverse, gen_smi, seen_smiles)
            is_lead_like = lead_like_enough(ll_filter, gen_mol)
            if is_new and is_lead_like:
                write_out(gen_mol, count, gen_smi, output)
                count += 1
    else:
        # fragmenting ---------------------------------------------------------
        mol_supplier = RobustSmilesMolSupplier(input_fn)
        for name, mol in mol_supplier:
            for i in range(nb_passes):
                fragments_smi, parent_name, dico = cut_some_bonds(mol, rng_seed)
                print("%s\t%s_p%d;%s" %
                      (fragments_smi, name, i, str(dico)), file=output)
            count += 1
    after = time.time()
    dt = after - before
    if assemble:
        print("generated %d molecules at %.2f mol/s" %
              (count, count / dt), file=sys.stderr)
    else:
        print("read %d molecules at %.2f mol/s" %
              (count, count / dt), file=sys.stderr)
    output.close()
