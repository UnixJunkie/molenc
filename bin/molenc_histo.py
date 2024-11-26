#!/usr/bin/env python3
#
# usage: ./histo.py FILE:str MIN:float MAX:float STEPS:int
#        0          1        2         3         4

import sys

input_fn = sys.argv[1]
min_val = float(sys.argv[2])
max_val = float(sys.argv[3])
assert(min_val < max_val)
num_steps = float(int(sys.argv[4]))

delta = (max_val - min_val) / num_steps

print('DEBUG: min:%g max:%g steps:%d delta:%g' %
      (min_val, max_val, num_steps, delta), file=sys.stderr)

# initialize the histogram
histo = {}
x = min_val
i = 0
while x < max_val:
    histo[i] = 0
    x += delta
    i += 1

# finalize the histogram
for line in open(input_fn).readlines():
    strip = line.strip()
    x = float(strip)
    assert(x >= min_val)
    assert(x <= max_val)
    hist_bin = int((x - min_val) / delta)
    histo[hist_bin] += 1

# print out the histogram
x = min_val
i = 0
while x < max_val:
    x_val = min_val + float(i) * delta
    y_val = histo[i]
    print('%f %d' % (x_val, y_val))
    x += delta
    i += 1
