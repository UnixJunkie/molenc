#!/bin/bash

#set -x # DEBUG

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
rm -f data/alcools.AD.curr
_build/default/src/AP_BBAD.exe -i data/alcools.smi -o data/alcools.AD.curr
diff data/alcools.AD.curr data/alcools.AD.ref

# 4) the BBAD of some molecules doesn't filter out any of those molecules
rm -f filtered.txt
_build/default/src/AP_BBAD.exe --bbad seq_AD.txt -i data/chembl1868_std.smi -o filtered.txt -np ${nprocs}
diff <(cat data/chembl1868_std.smi | wc -l) <(cat filtered.txt | wc -l)

# 5) the BBAD union for two sets of molecules should be the same as the AD for the union of the sets
rm -f head_AD.txt tail_AD.txt head_tail_AD_union.txt head_tail_AD.txt
_build/default/src/AP_BBAD.exe -i <(head data/chembl1868_std.smi) -o head_AD.txt -np ${nprocs}
_build/default/src/AP_BBAD.exe -i <(tail data/chembl1868_std.smi) -o tail_AD.txt -np ${nprocs}
_build/default/src/AP_BBAD.exe --bbad head_AD.txt,tail_AD.txt -o head_tail_AD_union.txt
_build/default/src/AP_BBAD.exe -i <(head data/chembl1868_std.smi; tail data/chembl1868_std.smi) \
                               -o head_tail_AD.txt
diff head_tail_AD_union.txt head_tail_AD.txt
