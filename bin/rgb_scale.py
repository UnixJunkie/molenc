#!/usr/bin/env python3

# pip3 install colour # to get the required library

from colour import Color

red = Color("red")
colors = list(red.range_to(Color("white"), 101))
for c in colors:
    (r, g, b) = c.get_rgb()
    print("%.2f %.2f %.2f" % (r, g, b))
