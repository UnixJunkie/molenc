#!/bin/bash

#set -x # DEBUG
set -u

if [ $# -eq 0 ]; then
    echo "usage:"
    echo "molenc.sh -i input.smi -o output.txt"
    echo "         [-d encoding.dix]: reuse existing feature dictionary"
    echo "         [-r i:j]: fingerprint radius (default=0:1)"
    echo "         [--pairs]: use atom pairs instead of Faulon's FP"
    echo "         [--seq]: sequential mode (disable parallelization)"
    echo "         [-v]: debug mode; keep temp files"
    echo "         [-n <int>]: max jobs in parallel"
    echo "         [-c <int>]: chunk size"
    echo "         [--no-std]: don't standardize input file molecules"
    echo "                     ONLY USE IF THEY HAVE ALREADY BEEN STANDARDIZED"
    exit 1
fi

input=""
output=""
dico=""
range="0:1" # default val
no_std=""
sequential=""
debug=""
nprocs="1"
csize="1"
pairs=""

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
        -n)
            nprocs="$2"
            shift
            shift
            ;;
        -c)
            csize="$2"
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
        --seq)
            sequential="TRUE"
            shift # past argument
            ;;
        --pairs)
            pairs="TRUE"
            shift # past argument
            ;;
        -v)
            debug="TRUE"
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

# enable parallelization or not
there_is_pardi=`which pardi`
if [ "$sequential" != "" ]; then
    there_is_pardi=""
fi

if [ "$no_std" == "" ]; then
    # tell user how to install standardiser if not here
    which standardiser 2>&1 > /dev/null || \
        echo 'molenc: ERROR: type: pip3 install chemo-standardizer'
    if [ $nprocs -gt 1 ] && [ "$there_is_pardi" != "" ]; then
        echo 'standardizing molecules in parallel...'
        pardi -p -n $nprocs -i $input -o $tmp_smi -c 100 -d l -ie '.smi' -oe '.smi' \
              -w 'standardiser -i %IN -o %OUT 2>/dev/null'
    else
        echo 'standardizing molecules...'
        (standardiser -i $input -o $tmp_smi 2>&1) > $std_log
    fi
else
    cp $input $tmp_smi
fi
if [ $nprocs -gt 1 ] && [ "$there_is_pardi" != "" ]; then
    echo 'typing atoms in parallel...'
    pardi -p -n $nprocs -i $tmp_smi -o $tmp_types -c 100 -d l -ie '.smi' \
          -w 'molenc_type_atoms.py -i %IN -o %OUT 2>/dev/null'
else
    echo 'typing atoms...'
    molenc_type_atoms.py -i $tmp_smi -o $tmp_types
fi

echo "encoding molecules..."
if [ "$pairs" != "" ]; then
    if [ "$dico" != "" ]; then
        molenc_ap -np $nprocs -cs $csize -i $tmp_types -o $output -id $dico
    else
        molenc_ap -np $nprocs -cs $csize -i $tmp_types -o $output -od $input'.dix'
    fi
else
    if [ "$dico" != "" ]; then
        # if dictionary is provided, parallelize encoding
        molenc_e -n $nprocs -i $tmp_types -r $range -o $tmp_enc -d $dico
        molenc_d -n $nprocs -i $tmp_enc -o $output -d $dico
    else
        molenc_e -i $tmp_types -r $range -o $tmp_enc
        molenc_d -i $tmp_enc -o $output
    fi
fi

# cleanup
if [ "$debug" == "" ]; then
    rm -f $std_log $tmp $tmp_smi $tmp_types $tmp_enc
fi
