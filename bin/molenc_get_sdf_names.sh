#!/bin/bash
#
# extract molecule names from properly formated SDF file
# usage: ./molenc_get_sdf_names.sh IN.sdf OUT.names

set -u

FN=$1
OUT=$2
head -1 ${FN} > ${OUT}
\egrep --no-group-separator -A1 '^\$\$\$\$$' ${FN} | \egrep -v '^\$\$\$\$$' >> ${OUT}
