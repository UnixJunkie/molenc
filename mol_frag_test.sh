#!/bin/bash

set -x

# clean
rm -f data/3.to_frag data/3_frags.txt data/3_frags.smi data/3_genmols.txt data/3_genmols.smi

# regen
./bin/molenc_frag.py -i data/3.smi -o data/3.to_frag
which mview && mview data/3.smi &
./molenc_frag -im data/3.to_frag -of data/3_frags.txt -s 1234
./bin/molenc_frag2smi.py -i data/3_frags.txt -o data/3_frags.smi
which mview && mview data/3_frags.smi &
./molenc_frag -if data/3_frags.txt -om data/3_genmols.txt -s 1234 -n 20
./bin/molenc_mol2smi.py -i data/3_genmols.txt -o data/3_genmols.smi
cut -f1 data/3_genmols.smi | sort -u > data/3_genmols_uniq.smi
which mview && mview data/3_genmols_uniq.smi &

# IN=data/all_kegg_drugs_20112019_std
# ./bin/molenc_frag.py -i $IN.smi -o $IN.to_frag
# ./molenc_frag -im $IN.to_frag -of $IN.frags -s 1234
# ./molenc_frag -if $IN.frags -om $IN.mols -s 1234 -n 50
