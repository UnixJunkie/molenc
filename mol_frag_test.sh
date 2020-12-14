#!/bin/bash

rm -f data/3.to_frag
./bin/molenc_frag.py -i data/3.smi -o data/3.to_frag
./molenc_frag -i data/3.to_frag -o data/3.frags -s 1234
