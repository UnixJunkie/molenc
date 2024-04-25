#!/bin/bash

# regression test for the UHD fingerprint
make

# cleanup any prior run
\rm -f data/ethanol.uhd data/ethanol.smi.dix

# run
_build/default/src/molenc_UHD.exe -f -i data/ethanol.smi -o data/ethanol.uhd

# check Vs refs
diff data/ethanol.uhd     data/ethanol.uhd.ref
diff data/ethanol.smi.dix data/ethanol.uhd.dix.ref
