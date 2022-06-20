#!/bin/bash

set -x # DEBUG

# check we parse correctly 10 3D molecules in a .sdf file
rm -f data/chembl30_10mols.txt.curr
_build/default/src/sdf_read.exe data/chembl30_10mols.sdf \
                                > data/chembl30_10mols.txt.curr
diff data/chembl30_10mols.txt.ref data/chembl30_10mols.txt.curr
