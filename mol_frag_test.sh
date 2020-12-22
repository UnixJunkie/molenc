#!/bin/bash

set -x

rm -f data/3.to_frag data/3.frags data/3.frags2 data/3.mols
./bin/molenc_frag.py -i data/3.smi -o data/3.to_frag
./molenc_frag -im data/3.to_frag -of data/3.frags -s 1234
./molenc_frag -im data/3.to_frag -of data/3.frags2 -s 1234 -p 2
./molenc_frag -if data/3.frags -om data/3.mols -s 1234 -n 3

# IN=data/all_kegg_drugs_20112019_std
# ./bin/molenc_frag.py -i $IN.smi -o $IN.to_frag
# ./bin/molenc_frag.py -i $IN.smi -o $IN.to_frag2 -p 2
# ./molenc_frag -im $IN.to_frag -of $IN.frags -s 1234
# ./molenc_frag -if $IN.frags -om $IN.mols -s 1234 -n 50
