#!/bin/bash

set -x #DEBUG

# clean
rm -f data/chembl_antivirals.frags.smi data/chembl_antivirals.genmols.txt
# fragment
./bin/molenc_smisur.py --seed 1234 \
                       -i data/chembl_antivirals.smi \
                       -o data/chembl_antivirals.frags.smi
# assemble
./bin/molenc_smisur.py --seed 1234 --assemble \
                       -i data/chembl_antivirals.frags.smi \
                       -o data/chembl_antivirals.genmols.txt
