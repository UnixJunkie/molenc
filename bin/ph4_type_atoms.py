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

# SMARTS feature families concatenated from
# /usr/share/RDKit/Data/BaseFeatures.fdef @ 18/10/2018
feat_to_smarts = {
    'A': '[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*),$([O;-;!$(*-N=O)]),$([o;+0])]),$([n;+0;!X3;!$([n;H1](cc)cc),$([$([N;H0]#[C&v4])]),$([N&v3;H0;$(Nc)])]),$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]',
    'a': '[$([a;r4,!R1&r3])]1:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:1,[$([a;r5,!R1&r4,!R1&r3])]1:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:1,[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:1,[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1,[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1',
    'D': '[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),$([$(n[n;H1]),$(nc[n;H1])])]),$([O,S;H1;+0])]',
    'H': '[R0;D2;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])],[D3,D4;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]',
    'h': '[N;D3;+](=O)[O-],[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1,[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1,[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1,[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1,[CH;!R](-[CH3])-[CH3],[C;!R](-[CH3])(-[CH3])-[CH3]',
    'N': '[C,S](=[O,S,P])-[O;H1,H0&-1]',
    'P': '[$([$([N;H2&+0][$([C;!$(C=*)])])]),$([$([N;H1&+0]([$([C;!$(C=*)])])[$([C;!$(C=*)])])]),$([$([N;H0&+0]([$([C;!$(C=*)])])([$([C;!$(C=*)])])[$([C;!$(C=*)])])]);!$(N[a])],NC(=N)N,c1ncnc1,[#7;+;!$([N+]-[O-])]',
    'Z': '[S;D1]-[#6],[#6]-C(=O)-C-[S;D1],[#6]-C(=O)-C-C-[S;D1],[#6]-C(=O)-N-[O;D1],[#6]-C(=O)-[O;D1],[#6]-P(=O)(-O)-[C,O,N]-[C,H]' }

acceptor_smarts =     feat_to_smarts['A']
aromatic_smarts =     feat_to_smarts['a']
donor_smarts =        feat_to_smarts['D']
hydro_smarts =        feat_to_smarts['H']
lumped_hydro_smarts = feat_to_smarts['h']
neg_smarts =          feat_to_smarts['N']
pos_smarts =          feat_to_smarts['P']
zn_smarts =           feat_to_smarts['Z']

acceptor_pattern =     Chem.MolFromSmarts(acceptor_smarts)
aromatic_pattern =     Chem.MolFromSmarts(aromatic_smarts)
donor_pattern =        Chem.MolFromSmarts(donor_smarts)
hydro_pattern =        Chem.MolFromSmarts(hydro_smarts)
lumped_hydro_pattern = Chem.MolFromSmarts(lumped_hydro_smarts)
neg_pattern =          Chem.MolFromSmarts(neg_smarts)
pos_pattern =          Chem.MolFromSmarts(pos_smarts)
zn_pattern =           Chem.MolFromSmarts(zn_smarts)

def add_to_set(tt, s):
    for i in tt:
        for j in i:
            s.add(j)
    return s

def get_ph4_feats(mol):
    acc = mol.GetSubstructMatches(acceptor_pattern)
    aro = mol.GetSubstructMatches(aromatic_pattern)
    don = mol.GetSubstructMatches(donor_pattern)
    hyd = mol.GetSubstructMatches(hydro_pattern)
    lhy = mol.GetSubstructMatches(lumped_hydro_pattern)
    neg = mol.GetSubstructMatches(neg_pattern)
    pos = mol.GetSubstructMatches(pos_pattern)
    zns = mol.GetSubstructMatches(zn_pattern)
    acc_indexes = add_to_set(acc, Set([]))
    aro_indexes = add_to_set(aro, Set([]))
    don_indexes = add_to_set(don, Set([]))
    hyd_indexes = add_to_set(hyd, Set([]))
    lhy_indexes = add_to_set(lhy, Set([]))
    neg_indexes = add_to_set(neg, Set([]))
    pos_indexes = add_to_set(pos, Set([]))
    zn__indexes = add_to_set(zns, Set([]))
    atom_index_to_features = {}
    for a in mol.GetAtoms():
        id = a.GetIdx()
        atom_index_to_features[id] = Set([])
    for i in acc_indexes:
        atom_index_to_features[i].add('A')
    for i in aro_indexes:
        atom_index_to_features[i].add('a')
    for i in don_indexes:
        atom_index_to_features[i].add('D')
    for i in hyd_indexes:
        atom_index_to_features[i].add('H')
    for i in lhy_indexes:
        atom_index_to_features[i].add('h')
    for i in neg_indexes:
        atom_index_to_features[i].add('N')
    for i in pos_indexes:
        atom_index_to_features[i].add('P')
    for i in zns_indexes:
        atom_index_to_features[i].add('Z')
    return atom_index_to_features

def get_mol_feats(mol):
    feats = get_ph4_feats(mol)
    for a in mol.GetAtoms():
      id = a.GetIdx()
      features = feats[i]
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

# factory.GetFeatureDefs()
# {'Acceptor.SingleAtomAcceptor': '[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*),$([O;-;!$(*-N=O)]),$([o;+0])]),$([n;+0;!X3;!$([n;H1](cc)cc),$([$([N;H0]#[C&v4])]),$([N&v3;H0;$(Nc)])]),$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]',
#  'Aromatic.Arom4': '[$([a;r4,!R1&r3])]1:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:1',
#  'Aromatic.Arom5': '[$([a;r5,!R1&r4,!R1&r3])]1:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:1',
#  'Aromatic.Arom6': '[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:1',
#  'Aromatic.Arom7': '[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1',
#  'Aromatic.Arom8': '[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1',
#  'Donor.SingleAtomDonor': '[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),$([$(n[n;H1]),$(nc[n;H1])])]),$([O,S;H1;+0])]',
#  'Hydrophobe.ChainTwoWayAttach': '[R0;D2;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]',
#  'Hydrophobe.ThreeWayAttach': '[D3,D4;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]',
#  'LumpedHydrophobe.Nitro2': '[N;D3;+](=O)[O-]',
#  'LumpedHydrophobe.RH3_3': '[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
#  'LumpedHydrophobe.RH4_4': '[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
#  'LumpedHydrophobe.RH5_5': '[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
#  'LumpedHydrophobe.RH6_6': '[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
#  'LumpedHydrophobe.iPropyl': '[CH;!R](-[CH3])-[CH3]',
#  'LumpedHydrophobe.tButyl': '[C;!R](-[CH3])(-[CH3])-[CH3]',
#  'NegIonizable.AcidicGroup': '[C,S](=[O,S,P])-[O;H1,H0&-1]',
#  'PosIonizable.BasicGroup': '[$([$([N;H2&+0][$([C;!$(C=*)])])]),$([$([N;H1&+0]([$([C;!$(C=*)])])[$([C;!$(C=*)])])]),$([$([N;H0&+0]([$([C;!$(C=*)])])([$([C;!$(C=*)])])[$([C;!$(C=*)])])]);!$(N[a])]',
#  'PosIonizable.Guanidine': 'NC(=N)N',
#  'PosIonizable.Imidazole': 'c1ncnc1',
#  'PosIonizable.PosN': '[#7;+;!$([N+]-[O-])]',
#  'ZnBinder.ZnBinder1': '[S;D1]-[#6]',
#  'ZnBinder.ZnBinder2': '[#6]-C(=O)-C-[S;D1]',
#  'ZnBinder.ZnBinder3': '[#6]-C(=O)-C-C-[S;D1]',
#  'ZnBinder.ZnBinder4': '[#6]-C(=O)-N-[O;D1]',
#  'ZnBinder.ZnBinder5': '[#6]-C(=O)-[O;D1]',
#  'ZnBinder.ZnBinder6': '[#6]-P(=O)(-O)-[C,O,N]-[C,H]'}

if __name__ == '__main__':
    before = time.time()
    fn = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    fdef_str = open_fn(fn).read()
    print(fdef_str)
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
