#!/usr/bin/env python3
#
# usage: ./histo.py FILE:str STEPS:int
#        0          1        2

import sys

if len(sys.argv) == 1:
    print("usage: molenc_histo.py INPUT_FILE NUM_STEPS", file=sys.stderr)
    exit(1)

input_fn = sys.argv[1]
num_steps = float(int(sys.argv[2]))

# read in values
floats = []
for line in open(input_fn).readlines():
    strip = line.strip()
    x = float(strip)
    floats.append(x)
assert(len(floats) > 0)

min_val = min(floats)
max_val = max(floats)
assert(min_val < max_val)
delta = (max_val - min_val) / num_steps
print('DEBUG: min:%g max:%g steps:%d delta:%g' %
      (min_val, max_val, num_steps, delta), file=sys.stderr)

# initialize the histogram
histo = {}
x = min_val
i = 0
while x <= max_val + delta:
    histo[i] = 0
    x += delta
    i += 1

# finalize the histogram
for x in floats:
    assert(x >= min_val)
    assert(x <= max_val)
    hist_bin = int((x - min_val) / delta)
    histo[hist_bin] += 1

# print out histogram
x = min_val
i = 0
while x <= max_val + delta:
    x_val = min_val + float(i) * delta
    y_val = histo[i]
    print('%f %d' % (x_val, y_val))
    x += delta
    i += 1
