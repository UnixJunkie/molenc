#!/usr/bin/env bash

set -u
set -x # DEBUG

# from SMILES to 3D protonated molecules with partial charges ready-to-dock

# FBR: check a .smi file was given as (only?) parameter

IN=$1
BASEDIR=`dirname $IN`
BASENAME=`basename $IN .smi`
OUT=${BASEDIR}'/'${BASENAME}

PROTONATED=${OUT}_taut74.smi
CONFORMER=${OUT}_taut74_1conf.sdf
CHARGED=${OUT}_taut74_1conf_mmff.mol2

# Default: fast mode ----------------------------------------------------------

# protonation state at physiological pH
obabel $IN -O $PROTONATED -p 7.4

# lowest energy conformer
~/usr/openeye/bin/omega2 -strictstereo false -maxconfs 1 \
                         -in $PROTONATED -out $CONFORMER

# assign partial charges
~/usr/openeye/bin/molcharge -method mmff -in $CONFORMER \
                            -out $CHARGED

# rm temp. files
rm -f $PROTONATED $CONFORMER

# user feedback
wc -l $IN
grep -c '@<TRIPOS>MOLECULE' $CHARGED

# if user really want: accurate model -----------------------------------------
# TODO
