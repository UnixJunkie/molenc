#!/usr/bin/env python

# usage: list_features.py molecules.sdf

from __future__ import print_function

import common
import os, sys, time
from sets import Set
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem

# we use a set of features so that there are no duplicated features for
# a given atom
def get_feats(dico, key):
    res = None
    try:
        res = dico[key]
    except KeyError:
        res = Set([])
    return res

def set_is_empty(s):
    s == Set([])

feature_family_to_char = { 'Donor': 'D',
                           'Acceptor': 'A',
                           'PosIonizable': 'P',
                           'NegIonizable': 'N',
                           'Aromatic': 'a',
                           'Hydrophobe': 'H',
                           'LumpedHydrophobe': 'h',
                           'ZnBinder': 'Z' }

def feat_to_char(feat):
    return feature_family_to_char[feat]

# FBR: check against factory.GetFeaturesDef()

acceptor = '[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*),$([O;-;!$(*-N=O)]),$([o;+0])]),$([n;+0;!X3;!$([n;H1](cc)cc),$([$([N;H0]#[C&v4])]),$([N&v3;H0;$(Nc)])]),$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]'
arom4 = '[$([a;r4,!R1&r3])]1:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:1'
arom5 = '[$([a;r5,!R1&r4,!R1&r3])]1:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:1'
arom6 = '[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:1'
arom7 = '[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1'
arom8 = '[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1'
donor = '[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),$([$(n[n;H1]),$(nc[n;H1])])]),$([O,S;H1;+0])]'
hydro1 = '[R0;D2;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]'
hydro2 = '[D3,D4;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]'
lhydro1 = '[N;D3;+](=O)[O-]'
lhydro2 = '[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1'
lhydro3 = '[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1'
lhydro4 = '[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1'
lhydro5 = '[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1'
lhydro6 = '[CH;!R](-[CH3])-[CH3]'
lhydro7 = '[C;!R](-[CH3])(-[CH3])-[CH3]'
neg = '[C,S](=[O,S,P])-[O;H1,H0&-1]'
pos1 = '[$([$([N;H2&+0][$([C;!$(C=*)])])]),$([$([N;H1&+0]([$([C;!$(C=*)])])[$([C;!$(C=*)])])]),$([$([N;H0&+0]([$([C;!$(C=*)])])([$([C;!$(C=*)])])[$([C;!$(C=*)])])]);!$(N[a])]'
pos2 = 'NC(=N)N'
pos3 = 'c1ncnc1'
pos4 = '[#7;+;!$([N+]-[O-])]'
Zn1 = '[S;D1]-[#6]'
Zn2 = '[#6]-C(=O)-C-[S;D1]'
Zn3 = '[#6]-C(=O)-C-C-[S;D1]'
Zn4 = '[#6]-C(=O)-N-[O;D1]'
Zn5 = '[#6]-C(=O)-[O;D1]'
Zn6 = '[#6]-P(=O)(-O)-[C,O,N]-[C,H]'

acc_pat = Chem.MolFromSmarts(acceptor)
arom4_pat = Chem.MolFromSmarts(arom4)
arom5_pat = Chem.MolFromSmarts(arom5)
arom6_pat = Chem.MolFromSmarts(arom6)
arom7_pat = Chem.MolFromSmarts(arom7)
arom8_pat = Chem.MolFromSmarts(arom8)
donor_pat = Chem.MolFromSmarts(donor)
hydro1_pat = Chem.MolFromSmarts(hydro1)
hydro2_pat = Chem.MolFromSmarts(hydro2)
lhydro1_pat = Chem.MolFromSmarts(lhydro1)
lhydro2_pat = Chem.MolFromSmarts(lhydro2)
lhydro3_pat = Chem.MolFromSmarts(lhydro3)
lhydro4_pat = Chem.MolFromSmarts(lhydro4)
lhydro5_pat = Chem.MolFromSmarts(lhydro5)
lhydro6_pat = Chem.MolFromSmarts(lhydro6)
lhydro7_pat = Chem.MolFromSmarts(lhydro7)
neg_pat = Chem.MolFromSmarts(neg)
pos1_pat = Chem.MolFromSmarts(pos1)
pos2_pat = Chem.MolFromSmarts(pos2)
pos3_pat = Chem.MolFromSmarts(pos3)
pos4_pat = Chem.MolFromSmarts(pos4)
Zn1_pat = Chem.MolFromSmarts(Zn1)
Zn2_pat = Chem.MolFromSmarts(Zn2)
Zn3_pat = Chem.MolFromSmarts(Zn3)
Zn4_pat = Chem.MolFromSmarts(Zn4)
Zn5_pat = Chem.MolFromSmarts(Zn5)
Zn6_pat = Chem.MolFromSmarts(Zn6)

def matching_indexes(mol, pat_str):
    res = []
    pat = mol.GetSubstructMatches(pat_str)
    for i in pat:
        for j in i:
            res.append(j)
    return res

def get_ph4_feats(mol):
    acc_match = matching_indexes(mol, acc_pat)
    arom4_match = matching_indexes(mol, arom4_pat)
    arom5_match = matching_indexes(mol, arom5_pat)
    arom6_match = matching_indexes(mol, arom6_pat)
    arom7_match = matching_indexes(mol, arom7_pat)
    arom8_match = matching_indexes(mol, arom8_pat)
    donor_match = matching_indexes(mol, donor_pat)
    hydro1_match = matching_indexes(mol, hydro1_pat)
    hydro2_match = matching_indexes(mol, hydro2_pat)
    lhydro1_match = matching_indexes(mol, lhydro1_pat)
    lhydro2_match = matching_indexes(mol, lhydro2_pat)
    lhydro3_match = matching_indexes(mol, lhydro3_pat)
    lhydro4_match = matching_indexes(mol, lhydro4_pat)
    lhydro5_match = matching_indexes(mol, lhydro5_pat)
    lhydro6_match = matching_indexes(mol, lhydro6_pat)
    lhydro7_match = matching_indexes(mol, lhydro7_pat)
    neg_match = matching_indexes(mol, neg_pat)
    pos1_match = matching_indexes(mol, pos1_pat)
    pos2_match = matching_indexes(mol, pos2_pat)
    pos3_match = matching_indexes(mol, pos3_pat)
    pos4_match = matching_indexes(mol, pos4_pat)
    zn1_match = matching_indexes(mol, Zn1_pat)
    zn2_match = matching_indexes(mol, Zn2_pat)
    zn3_match = matching_indexes(mol, Zn3_pat)
    zn4_match = matching_indexes(mol, Zn4_pat)
    zn5_match = matching_indexes(mol, Zn5_pat)
    zn6_match = matching_indexes(mol, Zn6_pat)
    atom_index_to_features = {}
    # create all needed sets, empty for the moment
    for a in mol.GetAtoms():
        id = a.GetIdx()
        atom_index_to_features[id] = Set([])
    for i in acc_match:
        prev = atom_index_to_features[i]
        prev.add('A')
        atom_index_to_features[i] = prev
    for arom in [arom4_match, arom5_match, arom6_match, arom7_match, arom8_match]:
        for i in arom:
            atom_index_to_features[i].add('a')
    for i in donor_match:
        atom_index_to_features[i].add('D')
    for hydro in [hydro1_match, hydro2_match]:
        for i in hydro:
            atom_index_to_features[i].add('H')
    for lhydro in [lhydro1_match, lhydro2_match, lhydro3_match, lhydro4_match, lhydro5_match, lhydro6_match, lhydro7_match]:
        for i in lhydro:
            atom_index_to_features[i].add('h')
    for i in neg_match:
        atom_index_to_features[i].add('N')
    for pos in [pos1_match, pos2_match, pos3_match, pos4_match]:
        for i in pos:
            atom_index_to_features[i].add('P')
    for zn in [zn1_match, zn2_match, zn3_match, zn4_match, zn5_match, zn6_match]:
        for i in zn:
            atom_index_to_features[i].add('Z')
    return atom_index_to_features

def get_mol_feats(mol):
    feats = get_ph4_feats(mol)
    for a in mol.GetAtoms():
      id = a.GetIdx()
      features = feats[id]
      if set_is_empty(features):
          print("%d" % id)
      else:
          str = " ".join(c for c in features)
          print("%d %s" % (id, str))

def ShowMolFeats(mol, factory):
  feat_map = {}
  molFeats = factory.GetFeaturesForMol(mol)
  for feat in molFeats:
    family = feat.GetFamily()
    # pos = feat.GetPos()
    # create atom_id to features set (some atoms don't have any ph4 feature)
    for id in feat.GetAtomIds():
        prev_feats = get_feats(feat_map, id)
        prev_feats.add(family)
        feat_map[id] = prev_feats
  for a in mol.GetAtoms():
      id = a.GetIdx()
      features = get_feats(feat_map, id)
      if set_is_empty(features):
          print("%d" % id)
      else:
          feat_chars = map(feat_to_char, features)
          str = " ".join(c for c in feat_chars)
          print("%d %s" % (id, str))

def open_fn(fn):
    res = None
    try:
        res = open(fn, 'r')
    except IOError:
        print('list_features.py: open_fn: could not open file %s' % fn,
              file = sys.stderr)
        sys.exit(1)
    return res

if __name__ == '__main__':
    before = time.time()
    fn = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    fdef_str = open_fn(fn).read()
    factory = AllChem.BuildFeatureFactoryFromString(fdef_str)
    mol_supplier = Chem.SDMolSupplier(sys.argv[1])
    count = 0
    for mol in mol_supplier:
        print("#mol %s" % mol.GetProp('_Name'))
        # ShowMolFeats(mol, factory)
        get_mol_feats(mol)
        common.print_bonds(mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
