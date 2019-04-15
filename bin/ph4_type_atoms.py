#!/usr/bin/env python3

# usage: list_features.py molecules.sdf

import common
import os, sys, time
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem

def open_fn(fn):
    res = None
    try:
        res = open(fn, 'r')
    except IOError:
        print('ph4_type_atoms.py: open_fn: could not open file %s' % fn,
              file = sys.stderr)
        sys.exit(1)
    return res

# create the features factory
fn = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
fdef_str = open_fn(fn).read()
factory = AllChem.BuildFeatureFactoryFromString(fdef_str)
feature_defs = factory.GetFeatureDefs()

## documentation
# feature_family_to_char = { 'Donor': 'D',
#                            'Acceptor': 'A',
#                            'PosIonizable': 'P',
#                            'NegIonizable': 'N',
#                            'Aromatic': 'a',
#                            'Hydrophobe': 'H',
#                            'LumpedHydrophobe': 'h',
#                            'ZnBinder': 'Z' }

# feature defs come from factory.GetFeaturesDef()
# I used all the keys found on 18/10/2018 DD/MM/YYYY
acceptor = feature_defs['Acceptor.SingleAtomAcceptor']
arom4 = feature_defs['Aromatic.Arom4']
arom5 = feature_defs['Aromatic.Arom5']
arom6 = feature_defs['Aromatic.Arom6']
arom7 = feature_defs['Aromatic.Arom7']
arom8 = feature_defs['Aromatic.Arom8']
donor = feature_defs['Donor.SingleAtomDonor']
hydro1 = feature_defs['Hydrophobe.ChainTwoWayAttach']
hydro2 = feature_defs['Hydrophobe.ThreeWayAttach']
lhydro1 = feature_defs['LumpedHydrophobe.Nitro2']
lhydro2 = feature_defs['LumpedHydrophobe.RH3_3']
lhydro3 = feature_defs['LumpedHydrophobe.RH4_4']
lhydro4 = feature_defs['LumpedHydrophobe.RH5_5']
lhydro5 = feature_defs['LumpedHydrophobe.RH6_6']
lhydro6 = feature_defs['LumpedHydrophobe.iPropyl']
lhydro7 = feature_defs['LumpedHydrophobe.tButyl']
neg = feature_defs['NegIonizable.AcidicGroup']
pos1 = feature_defs['PosIonizable.BasicGroup']
pos2 = feature_defs['PosIonizable.Guanidine']
pos3 = feature_defs['PosIonizable.Imidazole']
pos4 = feature_defs['PosIonizable.PosN']
Zn1 = feature_defs['ZnBinder.ZnBinder1']
Zn2 = feature_defs['ZnBinder.ZnBinder2']
Zn3 = feature_defs['ZnBinder.ZnBinder3']
Zn4 = feature_defs['ZnBinder.ZnBinder4']
Zn5 = feature_defs['ZnBinder.ZnBinder5']
Zn6 = feature_defs['ZnBinder.ZnBinder6']

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
    # we use a set of features so that there are no duplicated features for
    # a given atom
    for a in mol.GetAtoms():
        id = a.GetIdx()
        atom_index_to_features[id] = set([])
    for i in acc_match:
        atom_index_to_features[i].add('A')
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
      if len(features) == 0:
          print("%d _" % id) # '_' means no ph4 feat.
      else:
          l = list(features)
          l.sort() # canonicalize feats list for given atom
          str = " ".join(c for c in l)
          print("%d %s" % (id, str))

if __name__ == '__main__':
    before = time.time()
    mol_supplier = Chem.SDMolSupplier(sys.argv[1])
    count = 0
    for mol in mol_supplier:
        print("#atoms:%d %s" % (mol.GetNumAtoms(), mol.GetProp('_Name')))
        get_mol_feats(mol)
        common.print_bonds(mol)
        common.print_distance_matrix(mol)
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
