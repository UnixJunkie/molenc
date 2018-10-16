#!/usr/bin/env python

from __future__ import print_function

import sys, os, getopt
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem

def ShowMolFeats(mol, factory):
  molFeats = factory.GetFeaturesForMol(mol)
  for i, feat in enumerate(molFeats):
    family = feat.GetFamily()
    pos = feat.GetPos(-1)
    aidText = ' '.join([str(x + 1) for x in feat.GetAtomIds()])
    print('%s\t%.3f\t%.3f\t%.3f\t1.0\t# %s' %
          (family, pos.x, pos.y, pos.z, aidText))

if __name__=='__main__':
    fdef_fn = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    fdef = None
    try:
        fdef = open(fdef_fn, 'r').read()
    except IOError:
        logger.error('ERROR: Could not open fdef file %s' % options.fdef)
        sys.exit(1)
    factory = AllChem.BuildFeatureFactoryFromString(fdef)
    mol_supplier = Chem.SDMolSupplier(sys.argv[1])
    for m in mol_supplier:
        print("#Molecule: %s" % m.GetProp('_Name'))
        ShowMolFeats(m, factory)
    sys.exit(0)
