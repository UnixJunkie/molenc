#!/usr/bin/env python3
#
# compute histogram for data file (one value per line)

import argparse, sys

if __name__ == '__main__':
    default_steps = 50
    # <CLI parsing> -----------------------------------------------------------
    parser = argparse.ArgumentParser(
        description = "Output histogram values for input file")
    parser.add_argument("-i", metavar = "input.csv", dest = "input_fn",
                        help = "input file; one value per line")
    parser.add_argument("-o", metavar = "output.csv", dest = "output_fn",
                        help = "output file for gnuplot histeps")
    parser.add_argument("-n", metavar = "num_steps", dest = "num_steps",
                        help = "number of histogram steps (default=%d)" %
                        default_steps)
    parser.add_argument("-min", metavar = "min_val", dest = "min_val",
                        help = "minimum value (default=auto)",
                        default=sys.float_info.max, type=float)
    parser.add_argument("-max", metavar = "max_val", dest = "max_val",
                        help = "maximum value (default=auto)",
                        default=-sys.float_info.max, type=float)
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print("usage: molenc_histo.py -i INPUT.csv -o OUTPUT.csv",
              file=sys.stderr)
        exit(1)

    input_fn = args.input_fn
    output_fn = args.output_fn
    num_steps = float(int(args.num_steps))
    mini = args.min_val
    maxi = args.max_val
    # </CLI parsing> ----------------------------------------------------------

    with open(output_fn, 'w') as output:
        # read in values
        floats = []
        for line in open(input_fn).readlines():
            strip = line.strip()
            x = float(strip)
            floats.append(x)
        assert(len(floats) > 0)

        min_val = min(floats)
        max_val = max(floats)
        # if user provided bounds; we don't allow going over them
        if min_val < mini:
            print("FATAL: actual min < min_val: %f < %f" % (min_val, mini), file=sys.stderr)
            exit(1)
        if max_val > maxi:
            print("FATAL: actual max > max_val: %f > %f" % (max_val, maxi), file=sys.stderr)
            exit(1)
        # for consistent histograms, the user will specify bounds
        min_val = min(mini, min_val)
        max_val = max(maxi, max_val)

        assert(min_val < max_val)
        delta = (max_val - min_val) / num_steps
        # print('DEBUG: min:%g max:%g steps:%d delta:%g' %
        #       (min_val, max_val, num_steps, delta), file=output)

        # initialize histogram
        histo = {}
        x = min_val
        i = 0
        while x <= max_val + delta:
            histo[i] = 0
            x += delta
            i += 1

        # finalize histogram
        for x in floats:
            assert(x >= min_val)
            assert(x <= max_val)
            hist_bin = int((x - min_val) / delta)
            histo[hist_bin] += 1

        # print histogram
        x = min_val
        i = 0
        while x <= max_val + delta:
            x_val = min_val + float(i) * delta
            y_val = histo[i]
            print('%f %d' % (x_val, y_val), file=output)
            x += delta
            i += 1
