#!/usr/bin/env python

# usage: list_features.py molecules.sdf

from __future__ import print_function

import common, os, sys, time
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

# SMARTS features from /usr/share/RDKit/Data/BaseFeatures.fdef @ 18/10/2018

# # FBR: turn them into simpler definitions
# Acceptor
# Aromatic
# Donor
# Hydrophobe
# LumpedHydrophobe
# NegIonizable
# PosIonizable
# ZnBinder

{'Acceptor': '[$([O;H1;v2]),$([O;H0;v2;!$(O=N-*),$([O;-;!$(*-N=O)]),$([o;+0])]),$([n;+0;!X3;!$([n;H1](cc)cc),$([$([N;H0]#[C&v4])]),$([N&v3;H0;$(Nc)])]),$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]'
 'Aromatic': '[$([a;r4,!R1&r3])]1:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:[$([a;r4,!R1&r3])]:1,[$([a;r5,!R1&r4,!R1&r3])]1:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:[$([a;r5,!R1&r4,!R1&r3])]:1,[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r6,!R1&r5,!R1&r4,!R1&r3])]:1,[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1,[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]1:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:[$([a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3])]:1',
 'Donor': '[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),$([$(n[n;H1]),$(nc[n;H1])])]),$([O,S;H1;+0])]',
 'Hydrophobe.ChainTwoWayAttach': '[R0;D2;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]',
 'Hydrophobe.ThreeWayAttach': '[D3,D4;$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])]',
 'LumpedHydrophobe.Nitro2': '[N;D3;+](=O)[O-]',
 'LumpedHydrophobe.RH3_3': '[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
 'LumpedHydrophobe.RH4_4': '[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
 'LumpedHydrophobe.RH5_5': '[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
 'LumpedHydrophobe.RH6_6': '[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1[$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])][$([r6,!R1&r5,!R1&r4,!R1&r3]);$([R;$([c,s,S&H0&v2,Br,I,$([#6;+0;!$([#6;$([#6]~[#7,#8,#9])])])])])]1',
 'LumpedHydrophobe.iPropyl': '[CH;!R](-[CH3])-[CH3]',
 'LumpedHydrophobe.tButyl': '[C;!R](-[CH3])(-[CH3])-[CH3]',
 'NegIonizable.AcidicGroup': '[C,S](=[O,S,P])-[O;H1,H0&-1]',
 'PosIonizable.BasicGroup': '[$([$([N;H2&+0][$([C;!$(C=*)])])]),$([$([N;H1&+0]([$([C;!$(C=*)])])[$([C;!$(C=*)])])]),$([$([N;H0&+0]([$([C;!$(C=*)])])([$([C;!$(C=*)])])[$([C;!$(C=*)])])]);!$(N[a])]',
 'PosIonizable.Guanidine': 'NC(=N)N',
 'PosIonizable.Imidazole': 'c1ncnc1',
 'PosIonizable.PosN': '[#7;+;!$([N+]-[O-])]',
 'ZnBinder.ZnBinder1': '[S;D1]-[#6]',
 'ZnBinder.ZnBinder2': '[#6]-C(=O)-C-[S;D1]',
 'ZnBinder.ZnBinder3': '[#6]-C(=O)-C-C-[S;D1]',
 'ZnBinder.ZnBinder4': '[#6]-C(=O)-N-[O;D1]',
 'ZnBinder.ZnBinder5': '[#6]-C(=O)-[O;D1]',
 'ZnBinder.ZnBinder6': '[#6]-P(=O)(-O)-[C,O,N]-[C,H]'}

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
    print(fdef_str)
    factory = AllChem.BuildFeatureFactoryFromString(fdef_str)
    mol_supplier = Chem.SDMolSupplier(sys.argv[1])
    count = 0
    for mol in mol_supplier:
        print("#mol %s" % mol.GetProp('_Name'))
        ShowMolFeats(mol, factory)
        common.print_bonds(mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
