#!/usr/bin/env python

# usage: list_features.py molecules.sdf

from __future__ import print_function

import sys, os
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

# FBR: check this is the complete list in the rdkit configuration file
#      BaseFeatures.fdef
feature_family_to_char = { 'Donor': 'D',
                           'Acceptor': 'A',
                           'PosIonizable': 'P',
                           'NegIonizable': 'N',
                           'Aromatic': 'a',
                           'Hydrophobe': 'H',
                           'LumpedHydrophobe': 'h' }

def feat_to_char(feat):
    return feature_family_to_char[feat]

def ShowMolFeats(mol, factory):
  feat_map = {}
  molFeats = factory.GetFeaturesForMol(mol)
  for i, feat in enumerate(molFeats):
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
    fn = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    fdef_str = open_fn(fn).read()
    factory = AllChem.BuildFeatureFactoryFromString(fdef_str)
    mol_supplier = Chem.SDMolSupplier(sys.argv[1])
    for m in mol_supplier:
        print("#name %s" % m.GetProp('_Name'))
        print("#atoms")
        ShowMolFeats(m, factory)
        print("#bonds")
        for b in m.GetBonds():
            begin = b.GetBeginAtomIdx()
            end = b.GetEndAtomIdx()
            print("%d %d" % (begin, end))
    sys.exit(0)
