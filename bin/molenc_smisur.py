#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# The Smiling Surgeon: a doctor operating directly at the SMILES level

# FRAGMENTATION
#
# 1) OK cut some cuttable bonds of a molecule
# 2) OK try to save it as SMILES to see what we get (we get a mixture)
# 3) OK atom-type only atoms at the ends of bonds that were cut
# 4) OK use isotope numbers as keys in the (former opposite) atom type map
# 6) OK output this SMILES plus the (int -> atom_type) map as a string
#       plus the parent molecule name

# FRAGMENTS ASSEMBLY
#
# 1) read in all fragments
# 2) create the map ((start_atom_type, end_atom_type) -> fragment)
# 3) choose a seed fragment randomly
# 4) if remaining attachment points on it:
#        draw a compatible fragment
#        attach it
#        rec_call
#    else:
#        return current_mol

import argparse
import ast
import molenc_common as common
import random
import rdkit
import sys
import time

from molenc_common import RobustSmilesMolSupplier
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.AtomPairs import Pairs

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
    name = mol.GetProp("name")
    index_to_atom_type = dict_reverse_binding(atom_type_to_index)
    return (smi, name, index_to_atom_type)

# SMILES fragmentation
def cut_some_bonds(mol, seed):
    cuttable_bonds = [b.GetIdx() for b in common.find_cuttable_bonds(mol)]
    total_weight = Descriptors.MolWt(mol)
    # 150 Da: D. Rognan's suggested max fragment weight
    nb_frags = round(total_weight / 150)
    max_cuts = min(len(cuttable_bonds), nb_frags - 1)
    # print("mol %s; cut %d bonds" % (mol.GetProp("name"), max_cuts),
    #       file=sys.stderr)
    random.shuffle(cuttable_bonds)
    to_cut = cuttable_bonds[0:max_cuts]
    if len(to_cut) == 0:
        # molecule too small: not fragmented
        # still, we output it so that input and output SMILES files can be
        # visualized side-by-side
        smi = Chem.MolToSmiles(mol)
        name = mol.GetProp("name")
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
    res = []
    for (smi, frag_name, dico) in FragmentsSupplier(fn):
        res.append((smi, frag_name, dico))
    return res

def choose_one_random_fragment(all_frags):
    n = len(all_frags)
    i = random.randint(0, n - 1)
    return all_frags[i]

def count_uniq_fragment(all_frags):
    ss = set() # string set
    for (smi, _frag_name, _dico) in all_frags:
        ss.add(smi)
    return len(ss)

def index_fragments(frags):
    res = {}
    for smi, frag_name, dico in frags:
        frag_mol = Chem.MolFromSmiles(smi)
        # process each attachment point
        for a in frag_mol.GetAtoms():
            if a.GetAtomicNum() == 0: # wildcard atom
                isotope = a.GetIsotope()
                # print(isotope, dico) # debug
                dst_type = dico[isotope]
                neighbs = a.GetNeighbors()
                assert(len(neighbs) == 1)
                src_atom = neighbs[0]
                src_type = common.type_atom(src_atom)
                # record the fragment under key: (dst_type, src_type)
                # (i.e. ready to use by requiring fragment)
                key = (dst_type, src_type)
                new_val = (frag_mol, frag_name, dico)
                try:
                    previous_frags = res[key]
                    previous_frags.append(new_val)
                except KeyError:
                    res[key] = [new_val]
    return res

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules, or assemble molecular fragments")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--seed", dest = "seed", default = 1234,
                        type = int, help = "RNG seed")
    parser.add_argument("-n", dest = "nb_passes", default = 1,
                        type = int, help = "number of fragmentation passes")
    parser.add_argument("--assemble", dest = "assemble", action = 'store_true',
                        default = False,
                        help = "assemble instead of fragmenting")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    nb_passes = args.nb_passes
    assemble = args.assemble
    rng_seed = args.seed
    random.seed(rng_seed)
    output = open(args.output_fn, 'w')
    count = 0
    if assemble: # assembling fragments ---------------------------------------
        fragments = read_all_fragments(input_fn)
        nb_uniq = count_uniq_fragment(fragments)
        print('%d fragments indexed (uniq: %d)' % (len(fragments), nb_uniq))
        seed_frag = choose_one_random_fragment(fragments)
        print('seed_frag: %s' % str(seed_frag)) # debug
        index = index_fragments(fragments)
        print(len(index)) # debug
        # inspect the index (to debug)
        for k, v in index.items():
            print("k:%s -> %d frags" % (k, len(v)))
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
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
    output.close()
