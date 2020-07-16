#!/bin/bash

# csv format output by molenc_lizard.py to txt format expected by
# several molenc tools

awk -F',' \
    '{print $1","$2",[0:"$3";1:"$4";2:"$5";3:"$6";4:"$7";5:"$8";6:"$9"]"}' $1
