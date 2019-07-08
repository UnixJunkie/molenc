#!/bin/bash

#set -x # DEBUG

# check input_fn and output_fn are passed
if [ "$#" -ne 2 ] && [ "$#" -ne 3 ]; then
    echo "usage: molenc.sh input.smi output.txt [encoding.dix]"
    exit 1
fi

input=$1
output=$2
dico=""
if [ "$#" -eq 3 ]; then
    dico=$3
fi

std_log=$input'.std_log'
tmp=`mktemp`
tmp_smi=$tmp'_std.smi'
tmp_types=$tmp'_std.types'
tmp_enc=$tmp'_std.enc'

# tell user how to install standardiser if not here
which standardiser 2>&1 > /dev/null || \
    echo 'ERROR: type: pip3 install chemo-standardizer'

echo standardizing molecules...
(standardiser -i $input -o $tmp_smi 2>&1) > $std_log
echo typing atoms...
molenc_type_atoms.py $tmp_smi > $tmp_types
echo encoding molecules...
molenc_e -i $tmp_types -r 0:1 -o $tmp_enc
if [ "$dico" != "" ]; then
    molenc_d -i $tmp_enc -o $output -d $dico
else
    molenc_d -i $tmp_enc -o $output
fi

# cleanup
rm -f $std_log $tmp $tmp_smi $tmp_types $tmp_enc
