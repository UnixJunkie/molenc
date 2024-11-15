#!/usr/bin/env python3

# remove listed tags from a .sdf file
# usage: molenc_sdf_strip.py input.sdf "<TAG1>,<TAG2>,..." > stripped.sdf

import sys

sdf_fn = sys.argv[1]
tags = sys.argv[2] # coma-separated list of tags to remove

tags_list = tags.split(',')
tags_set = set(tags_list)

skip = False

def endswith_any(line, tset):
    for tag in tset:
        if line.endswith(tag):
            return True
    return False

for line in open(sdf_fn).readlines():
    stripped = line.strip()
    if endswith_any(stripped, tags_set):
        skip = True # the tag line itself
    elif skip:
        skip = False # the line after
    else:
        print(line, end='') # any other line
