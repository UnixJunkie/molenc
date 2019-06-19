#!/bin/bash

set -x

input=$1
output=$2

# FBR: check there are two params; print usage if not

tmp=`mktemp`
tmp_smi=$tmp'_std.smi'
tmp_types=$tmp'_std.types'
tmp_enc=$tmp'_std.enc'

(standardiser -i $input -o $tmp_smi 2>&1) > $input'.std_log'
molenc_type_atoms.py $tmp_smi > $tmp_types
molenc_e -i $tmp_types -r 0:1 -o $tmp_enc
molenc_d -i $tmp_enc -o $output
