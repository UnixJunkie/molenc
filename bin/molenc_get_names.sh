#!/bin/bash
#
# extract molecule names from properly formated SDF files to stdout
# usage: ./molenc_get_names.sh 1.sdf 2.sdf 3.sdf ...

set -u

function get_sdf_names () {
    (echo '$$$$'; cat $1) | \
        \egrep --no-group-separator -A1 '^\$\$\$\$$' | \egrep -v '^\$\$\$\$$'
}

for fn in "$@"; do
    case "$fn" in
        *.sdf.gz)
            zcat $fn | get_sdf_names /dev/stdin
            ;;
        *.sdf)
            get_sdf_names $fn
            ;;
        *) # unknown option
            echo "molenc_get_names.sh: unsupported file format: "$fn >> /dev/stderr
            exit 1
            ;;
    esac
done
