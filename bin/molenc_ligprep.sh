#!/usr/bin/env bash

set -u
#set -x # DEBUG

# from SMILES to 3D protonated molecules with partial charges ready-to-dock

# check input filename was given
if [ $# -eq 0 ]; then
    echo "usage:"
    echo "molenc_ligprep.sh molecules.smi"
    echo "         [--XXX]: XXX option"
    exit 1
fi

FAST_MODE=""
HIGH_ACCURACY=""
FREE_MODE="TRUE" # default

# parse CLI options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fast)
            FAST_MODE="TRUE"
            echo "fast mode"
            shift # past argument
            ;;
        --accurate)
            HIGH_ACCURACY="TRUE"
            echo "HA mode"
            shift # past argument
            ;;
        --free)
            FREE_MODE="TRUE"
            echo "free mode"
            shift # past argument
            ;;
        *) # unknown option
            echo "molenc_ligprep.sh: unknown option: "$1
            exit 1
            ;;
    esac
done

IN=$1
BASEDIR=`dirname $IN`
BASENAME=`basename $IN .smi`
OUT=${BASEDIR}'/'${BASENAME}

PROTONATED=${OUT}_taut74.smi
CONFORMER=${OUT}_taut74_1conf.sdf
CHARGED=${OUT}_taut74_1conf_mmff.mol2

# fast mode -------------------------------------------------------------------

if [ $FAST_MODE == "TRUE" ]; then
    # protonation state at physiological pH
    obabel $IN -O $PROTONATED -p 7.4
    # lowest energy conformer w/ OE omega
    ~/usr/openeye/arch/Ubuntu-18.04-x64/omega/omega2 \
        -strictstereo false -maxconfs 1 -in $PROTONATED -out $CONFORMER
    # assign partial charges
    ~/usr/openeye/arch/Ubuntu-18.04-x64/quacpac/molcharge \
        -method mmff -in $CONFORMER -out $CHARGED
fi

# high accuracy mode ----------------------------------------------------------

if [ $HIGH_ACCURACY == "TRUE" ]; then
    PROTONATED=${OUT}_taut74.sdf
    CONFORMER=${OUT}_taut74_1conf.sdf
    CHARGED=${OUT}_taut74_1conf_am1bcc.mol2
    # protonation state at physiological pH
    cp $IN $OUT.smiles # chemaxon cxcalc requires .smiles file extension...
    # -g --> skip erroneous molecules
    cxcalc -g majortautomer -H 7.4 -f sdf $OUT.smiles > $PROTONATED
    rm -f $OUT.smiles
    # lowest energy conformer w/ OE omega
    ~/usr/openeye/arch/Ubuntu-18.04-x64/omega/omega2 \
        -strictstereo false -maxconfs 1 -in $PROTONATED -out $CONFORMER
    # assign partial charges
    ~/usr/openeye/arch/Ubuntu-18.04-x64/quacpac/molcharge \
        -method am1bcc -in $CONFORMER -out $CHARGED
fi

# free / open-source mode -----------------------------------------------------

if [ $FREE_MODE == "TRUE" ]; then
    PROTONATED=${OUT}_taut74.smi
    CONFORMER=${OUT}_taut74_1conf.sdf
    CHARGED=${OUT}_taut74_1conf_mmff.mol2
    # protonation state at physiological pH and unsalt
    obabel $IN -O $PROTONATED -p 7.4 -r
    # -g --> skip erroneous molecules
    cxcalc -g majortautomer -H 7.4 -f sdf $OUT.smiles > $PROTONATED
    rm -f $OUT.smiles
    # lowest energy conformer
    obabel $PROTONATED -O $CONFORMER --gen3D
    # assign partial charges
    obabel $CONFORMER -O $CHARGED --partialcharge mmff94
fi

# cleanup ---------------------------------------------------------------------

# rm temp. files
rm -f $PROTONATED $CONFORMER

# user feedback
wc -l $IN
grep -c '@<TRIPOS>MOLECULE' $CHARGED
