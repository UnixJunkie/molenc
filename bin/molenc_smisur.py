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

def read_flat_fragments(flat_frags_fn):
    res = []
    with open(flat_frags_fn) as f:
        for line in f:
            smi, name_dict = line.strip().split("\t") # enforce TAB-separated
            name, dict_str = name_dict.split(";")
            dico = ast.literal_eval(dict_str)
            if len(dico) > 0: # molecule _was_ fragmented (not too small)
                res.append((smi, name, dico))
    return res

def read_all_fragments(in_frags_fn, flat_frags_out_fn):
    res = []
    with open(flat_frags_out_fn, 'w') as out:
        for (smi, name, dico) in FragmentsSupplier(in_frags_fn):
            res.append((smi, name, dico))
            out.write('%s\t%s;%s\n' % (smi, name, str(dico)))
    return res

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
                a.SetProp("src_typ", src_typ)
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
def bind_molecules(m1, m2, dst_idx, idx2):
    # print('m1: %s' % Chem.MolToSmiles(m1)) #debug
    # print('m2: %s' % Chem.MolToSmiles(m2)) #debug
    n1 = m1.GetNumAtoms()
    n2 = m2.GetNumAtoms()
    m = n1 + n2
    rw_mol = Chem.RWMol(Chem.CombineMols(m1, m2))
    assert(rw_mol.GetNumAtoms() == m)
    name1 = get_name(m1)
    name2 = get_name(m2)
    new_name = '%s|%s' % (name1, name2)
    set_name(rw_mol, new_name)
    ai = rw_mol.GetAtomWithIdx(dst_idx)
    assert(ai.GetAtomicNum() == 0) # attachment point
    dst_typ = ai.GetProp("dst_typ")
    src_idx = get_src_atom_idx(ai)
    src_typ = ai.GetProp("src_typ")
    dst_idx2 = n1 + idx2
    aj = rw_mol.GetAtomWithIdx(dst_idx2)
    assert(aj.GetAtomicNum() == 0) # attachment point
    dst_typ2 = aj.GetProp("dst_typ")
    src_idx2 = get_src_atom_idx(aj)
    src_typ2 = aj.GetProp("src_typ")
    if (dst_typ == src_typ2 and
        src_typ == dst_typ2):
        # attach. points are compatible
        rw_mol.AddBond(src_idx, src_idx2, Chem.rdchem.BondType.SINGLE)
        # remove former attachment points
        rw_mol.RemoveAtom(dst_idx2) # to not shift lower atom indexes
        rw_mol.RemoveAtom(dst_idx)
        return rw_mol
    else:
        # attach. points are incompatible !!!
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
        src_typ = dst_a.GetProp("src_typ")
        # draw compatible fragment
        key = (src_typ, dst_typ) # current to wanted direction
        # print('want key: %s' % str(key)) # debug
        possible_compat_frags = frags_index[key]
        compat_frag = random_choose_one(possible_compat_frags)
        (frag_mol2, dst_idx2) = compat_frag
        dst_a2 = frag_mol2.GetAtomWithIdx(dst_idx2)
        dst_typ2 = dst_a2.GetProp("dst_typ")
        src_typ2 = dst_a2.GetProp("src_typ")
        # check fragments compatibility
        assert(src_typ == dst_typ2)
        assert(dst_typ == src_typ2)
        # connect them
        new_mol = bind_molecules(frag_seed_mol, frag_mol2, dst_idx, dst_idx2)
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
        if not (gen_smi in seen_smiles):
            seen_smiles.add(gen_smi)
            return True
        else:
            return False

def lead_like_enough(ll_filter, mol):
    return ((not ll_filter) or lead_like_filter(mol))

# Tran-Nguyen, V. K., Jacquemard, C., & Rognan, D. (2020).
# LIT-PCBA: An unbiased data set for machine learning and virtual screening.
# Journal of chemical information and modeling, 60(9), 4263-4273.
def drug_like_filter(mol):
    MolW = Descriptors.MolWt(mol)
    if MolW <= 150 or MolW >= 800: # 150 < MolW < 800 Da
        return False
    cLogP = Descriptors.MolLogP(mol)
    if cLogP <= -3.0 or cLogP >= 5.0: # −3.0 < AlogP < 5.0
        return False
    RotB = Descriptors.NumRotatableBonds(mol)
    if RotB >= 15: # RotB < 15
        return False
    HBA = Descriptors.NumHAcceptors(mol)
    if HBA >= 10: # HBA < 10
        return False
    HBD = Descriptors.NumHDonors(mol)
    if HBD >= 10: # HBD < 10
        return False
    FC = Chem.rdmolops.GetFormalCharge(mol)
    if FC <= -2 or FC >= 2: # −2.0 < FC < 2.0
        return False
    return True # Still here? Drug-like then!

def drug_like_enough(dl_filter, mol):
    return ((not dl_filter) or drug_like_filter(mol))

# Lisurek, M., Rupp, B., Wichard, J., Neuenschwander, M., von Kries, J. P.,
# Frank, R., ... & Kühne, R. (2010).
# Design of chemical libraries with potentially bioactive molecules applying
# a maximum common substructure concept. Molecular diversity, 14(2), 401-408.
# SMARTS patterns kindly provided by Michael Lisurek
pat1 = Chem.MolFromSmarts('[C,c]S(=O)(=O)[F,Cl,Br,I]') # sulfonylhalide
pat2 = Chem.MolFromSmarts('[C,c]S(=O)(=O)O[CX4]') # sulfone_ester
pat3 = Chem.MolFromSmarts('C(=O)[F,Cl,Br,I]') # acylhalide
pat4 = Chem.MolFromSmarts('O=COC=O') # acidanhydride
pat5 = Chem.MolFromSmarts('c1([F,Cl,Br,I])ncccn1') # 2-halo_pyrimidine
pat6 = Chem.MolFromSmarts('[H]C=O') # aldehyde
pat7 = Chem.MolFromSmarts('C(=O)C(=O)') # 1,2-dicarbonyl
pat8 = Chem.MolFromSmarts('C1OC1') # epoxide
pat9 = Chem.MolFromSmarts('C1NC1') # aziridine
pat10 = Chem.MolFromSmarts('C(=O)S') # thioester
pat11 = Chem.MolFromSmarts('[#7]!@[#7]') # hydrazine
pat12 = Chem.MolFromSmarts('C=[CH2]') # ethenes
pat13 = Chem.MolFromSmarts('[H,*,!N][N;!R]=[C;!R]([*,H])[*,H]') # imine
pat14 = Chem.MolFromSmarts('[CX4]I') # alkyl_iodide
pat15 = Chem.MolFromSmarts('[Se]') # selenide
pat16 = Chem.MolFromSmarts('O-O') # peroxide
pat17 = Chem.MolFromSmarts('[NX3]!@[OX2]') # hetero-hetero_single_bond
pat18 = Chem.MolFromSmarts('[NX3]!@[NX3]') # hetero-hetero_single_bond
pat19 = Chem.MolFromSmarts('[NX3]!@[SX2]') # hetero-hetero_single_bond
pat20 = Chem.MolFromSmarts('[SX2]!@[SX2]') # hetero-hetero_single_bond
pat21 = Chem.MolFromSmarts('[SX2]!@[OX2]') # hetero-hetero_single_bond

def stable_filter(mol):
    if (mol.HasSubstructMatch(pat1) or
        mol.HasSubstructMatch(pat2) or
        mol.HasSubstructMatch(pat3) or
        mol.HasSubstructMatch(pat4) or
        mol.HasSubstructMatch(pat5) or
        mol.HasSubstructMatch(pat6) or
        mol.HasSubstructMatch(pat7) or
        mol.HasSubstructMatch(pat8) or
        mol.HasSubstructMatch(pat9) or
        mol.HasSubstructMatch(pat10) or
        mol.HasSubstructMatch(pat11) or
        mol.HasSubstructMatch(pat12) or
        mol.HasSubstructMatch(pat13) or
        mol.HasSubstructMatch(pat14) or
        mol.HasSubstructMatch(pat15) or
        mol.HasSubstructMatch(pat16) or
        mol.HasSubstructMatch(pat17) or
        mol.HasSubstructMatch(pat18) or
        mol.HasSubstructMatch(pat19) or
        mol.HasSubstructMatch(pat20) or
        mol.HasSubstructMatch(pat21)):
        return False
    else:
        return True

def stable_enough(s_filter, mol):
    return ((not s_filter) or stable_filter(mol))

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "fragment molecules, or assemble molecular fragments")
    parser.add_argument("-i", metavar = "input.smi", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.smi", dest = "output_fn",
                        help = "output file")
    parser.add_argument("--frag-dump", metavar = "frags_dump.smi",
                        dest = "flat_frags_out_fn",
                        help = "flat fragments dump file",
                        default = None)
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
    parser.add_argument("--drug-like", dest = "dl_filter", default = False,
                        action = "store_true",
                        help = "only generate drug-like molecules")
    parser.add_argument("--stable", dest = "s_filter", default = False,
                        action = "store_true",
                        help = "only generate stable molecules")
    # parse CLI
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    flat_frags_out_fn = args.flat_frags_out_fn
    nb_passes = args.nb_passes
    nmols = args.nmols
    assemble = nmols > 0
    rng_seed = args.seed
    diverse = args.diverse
    ll_filter = args.ll_filter
    dl_filter = args.dl_filter
    s_filter = args.s_filter
    seen_smiles = set()
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    output = open(args.output_fn, 'w')
    drug_filter_fails = 0
    lead_filter_fails = 0
    stable_filter_fails = 0
    not_new_fails = 0
    count = 0
    if assemble: # assembling fragments ---------------------------------------
        smi_fragments = read_all_fragments(input_fn, "/dev/null")
        nb_uniq = count_uniq_fragment(smi_fragments)
        print('read %d fragments (uniq: %d)' % (len(smi_fragments), nb_uniq))
        index = index_fragments(smi_fragments)
        fragments = extract_fragments(index)
        print('%d fragment keys' % len(index))
        # # inspect the index (to debug)
        # for k, v in index.items():
        #     print("k:%s -> %d frags" % (k, len(v)))
        while count < nmols:
            # FBR: parallelize here
            seed_frag = random_choose_one(fragments)
            # print('seed_frag: %s' % get_name(seed_frag)) # debug
            gen_mol = grow_fragment(seed_frag, index)
            gen_smi = Chem.MolToSmiles(gen_mol)
            is_new = new_enough(diverse, gen_smi, seen_smiles)
            if not is_new:
                not_new_fails += 1
            is_lead_like = lead_like_enough(ll_filter, gen_mol)
            if not is_lead_like:
                lead_filter_fails += 1
            is_drug_like = drug_like_enough(dl_filter, gen_mol)
            if not is_drug_like:
                drug_filter_fails += 1
            is_stable = stable_enough(s_filter, gen_mol)
            if not is_stable:
                stable_filter_fails += 1
            if is_new and is_lead_like and is_drug_like and is_stable:
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
        # log failures
        print("Fails: drug: %d lead: %d stable: %d new: %d" %
              (drug_filter_fails,
               lead_filter_fails,
               stable_filter_fails,
               not_new_fails))
    else:
        print("read %d molecules at %.2f mol/s" %
              (count, count / dt), file=sys.stderr)
    output.close()
    if flat_frags_out_fn != None:
        # read all fragments back and store them in a "flat" format
        _smi_fragments = read_all_fragments(args.output_fn, flat_frags_out_fn)
