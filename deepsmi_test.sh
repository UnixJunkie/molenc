#!/bin/bash

#set -x

head ~/src/FMGO/data/TCM_20k.smi > input.smi
dos2unix input.smi
./bin/molenc_deepsmi.py --no-rings    -i input.smi   -o output.dsmi
./bin/molenc_deepsmi.py --no-rings -d -i output.dsmi -o output.smi
diff input.smi output.smi
rm -f output.{dsmi,smi}

./bin/molenc_deepsmi.py --no-branches    -i input.smi   -o output.dsmi
./bin/molenc_deepsmi.py --no-branches -d -i output.dsmi -o output.smi
diff input.smi output.smi
rm -f output.{dsmi,smi}

./bin/molenc_deepsmi.py --no-branches --no-rings    -i input.smi   -o output.dsmi
./bin/molenc_deepsmi.py --no-branches --no-rings -d -i output.dsmi -o output.smi
diff input.smi output.smi
rm -f output.{dsmi,smi}
