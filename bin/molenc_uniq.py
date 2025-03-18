#!/usr/bin/env python3
#
# only print on stdout a line if its key was never seen before

import sys

input_fn = sys.argv[1]
sep = sys.argv[2]
# user-provided field number, as in awk, start at 1
field = int(sys.argv[3]) - 1

seen = {}

for line in open(input_fn).readlines():
    strip = line.strip()
    toks = strip.split(sep)
    key = toks[field]
    already_seen = seen.get(key, False)
    if not already_seen:
        print(line)
        seen[key] = True
