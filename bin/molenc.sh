#!/bin/bash

#set -x # DEBUG
set -u

if [ $# -eq 0 ]; then
    echo "usage:"
    echo "molenc.sh -i input.smi -o output.txt"
    echo "         [-d encoding.dix]; reuse existing dictionary"
    echo "         [-r i:j]; fingerprint radius (default=0:1)"
    echo "         [--no-std]; don't standardize input file molecules"
    echo "                     only use if they have already been standardized"
    exit 1
fi

input=""
output=""
dico=""
range="0:1" # default val
no_std=""

# parse CLI options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i)
            input="$2"
            shift # skip option specifier
            shift # skip option value
            ;;
        -o)
            output="$2"
            shift
            shift
            ;;
        -d)
            dico="$2"
            shift
            shift
            ;;
        -r)
            range="$2"
            shift
            shift
            ;;
        --no-std)
            no_std="TRUE"
            shift # past argument
            ;;
        *) # unknown option
            echo "molenc: unknown option: "$1
            exit 1
            ;;
    esac
done

std_log=$input'.std_log'
tmp=`mktemp`
tmp_smi=$tmp'_std.smi'
tmp_types=$tmp'_std.types'
tmp_enc=$tmp'_std.enc'

# FBR: parallelize standardization if pardi is detected
# FBR: parallelize atom typing if pardi is detected

if [ "$no_std" == "" ]; then
    # tell user how to install standardiser if not here
    which standardiser 2>&1 > /dev/null || \
        echo 'ERROR: type: pip3 install chemo-standardizer'
    echo standardizing molecules...
    (standardiser -i $input -o $tmp_smi 2>&1) > $std_log
else
    cp $input $tmp_smi
fi
echo typing atoms...
molenc_type_atoms.py $tmp_smi > $tmp_types
echo encoding molecules...
molenc_e -i $tmp_types -r $range -o $tmp_enc
if [ "$dico" != "" ]; then
    molenc_d -i $tmp_enc -o $output -d $dico
else
    molenc_d -i $tmp_enc -o $output
fi

# cleanup
rm -f $std_log $tmp $tmp_smi $tmp_types $tmp_enc
