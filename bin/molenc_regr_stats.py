#!/usr/bin/env python3
#
# output R2 and RMSE regression statistics from a file containing pairs
# of float values (one blank-separated pair per line)

import sklearn, sys

from sklearn.metrics import r2_score, root_mean_squared_error

def read_pair(line):
    strip = line.strip()
    tokens = line.split()
    x = float(tokens[0])
    y = float(tokens[1])
    xy = (x, y)
    return xy

if __name__ == '__main__':
    input_fn = sys.argv[1]
    xs = []
    ys = []
    with open(input_fn) as input:
        for line in input.readlines():
            x, y = read_pair(line)
            xs.append(x)
            ys.append(y)
    r2 = r2_score(xs, ys)
    rmse = root_mean_squared_error(xs, ys)
    print('R2=%.3f RMSE=%.3f' % (r2, rmse))
