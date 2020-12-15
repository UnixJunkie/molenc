#!/bin/bash

rm -f data/3.to_frag data/3.frags data/3.mols
./bin/molenc_frag.py -i data/3.smi -o data/3.to_frag
./molenc_frag -im data/3.to_frag -of data/3.frags -s 1234
./molenc_frag -if data/3.frags -om data/3.mols -s 1234 -n 3
