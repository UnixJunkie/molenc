#!/bin/bash

#set -x # DEBUG

# check params
if [ "$#" -ne 3 ]; then
    #            0         1         2          3
    echo "usage: molenc.sh input.smi output.txt features.dix"
    exit 1
fi

input=$1
output=$2
dico=$3

std_log=$input'.std_log'
tmp=`mktemp`
tmp_smi=$tmp'_std.smi'
tmp_scan=$tmp'_scan.smi'
tmp_types=$tmp'_std.types'
tmp_enc=$tmp'_std.enc'

# tell user how to install standardiser if not here
which standardiser 2>&1 > /dev/null || \
    echo 'ERROR: type: pip3 install chemo-standardizer'

echo standardizing molecules...
(standardiser -i $input -o $tmp_smi 2>&1) > $std_log
echo wildcard scan...
molenc_scan.py $tmp_smi > $tmp_scan
echo typing atoms...
molenc_type_atoms.py $tmp_scan > $tmp_types
echo encoding molecules...
molenc_e -i $tmp_types -r 0:1 -o $tmp_enc
molenc_d -i $tmp_enc -o $output -d $dico

# cleanup
rm -f $std_log $tmp $tmp_smi $tmp_scan $tmp_types $tmp_enc
