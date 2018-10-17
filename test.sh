#!/bin/bash

set -x # DEBUG

# encoding an SDF or a SMILES file is the same
diff <(./bin/type_atoms.py data/caffeine.sdf) \
     <(./bin/type_atoms.py data/caffeine.smi)

# ph4 features are the same than the ones extracted by ShowFeats.py
# (that were checked by hand and stored in a reference file)
diff <(./bin/ph4_type_atoms.py data/caffeine.sdf) data/caffeine_feats.ref
