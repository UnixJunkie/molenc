#!/usr/bin/env bash

set -u
#set -x # DEBUG

# from SMILES to 3D protonated molecules with partial charges ready-to-dock

# check input filename was given
if [ $# -eq 0 ]; then
    echo "usage:"
    echo "molenc_ligprep.sh molecules.smi"
    echo "         [--fast]: protonate:obabel; 3D:omega; molcharge:MMFF94"
    echo "         [--HA]: protonate:cxcalc; 3D:omega; molcharge:AM1BCC"
    echo "         [--free]: protonate:obabel; 3D:obabel; obabel:MMFF94"
    exit 1
fi

FAST_MODE=""
HIGH_ACCURACY=""
FREE_MODE="TRUE" # default

IN=$1
BASEDIR=`dirname $IN`
BASENAME=`basename $IN .smi`
OUT=${BASEDIR}'/'${BASENAME}

RHEL_OE_BASE=~/usr/openeye/arch/redhat-RHEL7-x64
UBUN_OE_BASE=~/usr/openeye/arch/Ubuntu-18.04-x64

# parse CLI options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fast)
            FAST_MODE="TRUE"
            HIGH_ACCURACY=""
            FREE_MODE=""
            echo "selected fast mode"
            shift # past argument
            ;;
        --HA)
            HIGH_ACCURACY="TRUE"
            FAST_MODE=""
            FREE_MODE=""
            echo "selected HA mode"
            shift # past argument
            ;;
        --free)
            FREE_MODE="TRUE"
            FAST_MODE=""
            HIGH_ACCURACY=""
            echo "selected free mode"
            shift # past argument
            ;;
        *) # unknown option or input file
            echo "input file: "$1
            shift # past argument
            ;;
    esac
done

PROTONATED=${OUT}_taut74.smi
CONFORMER=${OUT}_taut74_1conf.sdf
CHARGED=${OUT}_taut74_1conf_mmff.mol2

# fast mode -------------------------------------------------------------------

if [ "$FAST_MODE" == "TRUE" ]; then
    echo "running fast mode"
    # protonation state at physiological pH
    obabel $IN -O $PROTONATED -p 7.4
    # lowest energy conformer w/ OE omega
    ${RHEL_OE_BASE}/omega/omega2 \
        -strictstereo false -maxconfs 1 -in $PROTONATED -out $CONFORMER
    # assign partial charges
    ${RHEL_OE_BASE}/quacpac/molcharge \
        -method mmff -in $CONFORMER -out $CHARGED
fi

# high accuracy mode ----------------------------------------------------------

if [ "$HIGH_ACCURACY" == "TRUE" ]; then
    echo "running HA mode"
    PROTONATED=${OUT}_taut74.sdf
    CONFORMER=${OUT}_taut74_1conf.sdf
    CHARGED=${OUT}_taut74_1conf_am1bcc.mol2
    # protonation state at physiological pH
    cp $IN $OUT.smiles # chemaxon cxcalc requires .smiles file extension...
    # -g --> skip erroneous molecules
    cxcalc -g majortautomer -H 7.4 -f sdf $OUT.smiles > $PROTONATED
    rm -f $OUT.smiles
    # lowest energy conformer w/ OE omega
    ${RHEL_OE_BASE}/omega/omega2 \
        -strictstereo false -maxconfs 1 -in $PROTONATED -out $CONFORMER
    # assign partial charges
    ${RHEL_OE_BASE}/quacpac/molcharge \
        -method am1bcc -in $CONFORMER -out $CHARGED
fi

# free / open-source mode -----------------------------------------------------

if [ "$FREE_MODE" == "TRUE" ]; then
    echo "running free mode"
    PROTONATED=${OUT}_taut74.smi
    CONFORMER=${OUT}_taut74_1conf.sdf
    CHARGED=${OUT}_taut74_1conf_mmff.mol2
    # protonation state at physiological pH and unsalt
    obabel $IN -O $PROTONATED -p 7.4 -r
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
