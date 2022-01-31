#!/bin/bash

# check some properties of the AP-BBAD

# 1) the BBAD of a single molecule is the encoded single molecule
rm -f caffeine_AP_BBAD.txt
_build/default/src/AP_BBAD.exe -i data/caffeine.smi -o caffeine_AP_BBAD.txt
awk -v sum=0 -F' ' '{sum += $2}END{if(sum == 105){print "|features| OK"}}' caffeine_AP_BBAD.txt

# 2) the BBAD computed in parallel is the same as the sequential one
rm -f seq_AD.txt par_AD.txt
_build/default/src/AP_BBAD.exe -i data/chembl1868_std.smi -o seq_AD.txt -np 1
nprocs=`getconf _NPROCESSORS_ONLN`
_build/default/src/AP_BBAD.exe -i data/chembl1868_std.smi -o par_AD.txt -np ${nprocs}
diff seq_AD.txt par_AD.txt

# 3) compute a simple BBAD by hand; check this is the one we obtain

# 4) the BBAD of some molecules doesn't filter out any of those molecules

# 5) the BBAD union for two single molecules should be the same as 3)
