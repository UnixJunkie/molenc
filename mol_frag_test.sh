#!/bin/bash

rm -f data/3.to_frag
./bin/molenc_frag.py -i data/3.smi -o data/3.to_frag
./molenc_frag -im data/3.to_frag -of data/3.frags -s 1234
