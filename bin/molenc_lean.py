#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# Short: tool to automatically find the necessary/minimal training set size

# Long: find the size of the random partition from the training set that is
# enough to train a model (can help avoiding overfitting, since we might not
# need to use the whole training set to tune the model).
# However, note that for production you will still need to train your
# model on the full training set.
# In effect, we check that the distribution of the variable to model does not
# significantly change by adding more training samples.
# Usually, with a smaller dataset, you can get a model faster and hence
# accelerate your computational experiments.
#
# Bibliography:
# =============
# Domingos, P. (2012).
# "A few useful things to know about machine learning."
# Communications of the ACM, 55(10), 78-87.

import argparse, random, scipy, sys
import numpy as np
from scipy import stats

stderr = sys.stderr

def lines_of_file(fn):
    with open(fn) as f:
        return f.readlines()

def get_field(sep_char, field, line):
    tokens = line.split(sep = sep_char)
    value_str = tokens[field - 1]
    return float(value_str)

def get_all_values(sep_char, field, lines):
    res = []
    for l in lines:
        value = get_field(sep_char, field, l)
        res.append(value)
    return res

# ANSI terminal colors for UNIX
black   = "\033[30m"
red     = "\033[31m"
green   = "\033[32m"
yellow  = "\033[33m"
blue    = "\033[34m"
magenta = "\033[35m"
cyan    = "\033[36m"
white   = "\033[37m"
color_reset = "\033[39m"

def dichotomic_search(all_values, good_vals, low, curr, high):
    smaller_sample = np.random.choice(all_values, size = curr, replace = False)
    ks, p_val = stats.ks_2samp(smaller_sample, all_values,
                               alternative='two-sided', mode='auto')
    if p_val > 0.95:
        # recurse lower
        # print("%sbig_enough: %d" % (green, curr))
        good_vals.append(curr)
        curr_next = round((low + curr) / 2)
        if curr_next != curr:
            dichotomic_search(all_values, good_vals, low, curr_next, curr)
        else:
            return
    else:
        # recurse bigger
        # print("%stoo small: %d" % (white, curr))
        curr_next = round((curr + high) / 2)
        if curr_next != curr:
            dichotomic_search(all_values, good_vals, curr, curr_next, high)
        else:
            return

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "automatically determine minimal training set size")
    parser.add_argument("-i", metavar = "INPUT_FILE", dest = "input_fn",
                        help = "line-oriented file containing the \
                        target variable (training-set)", type = str)
    parser.add_argument("-d", metavar = "DELIM_CHAR", dest = "sep_char",
                        help = "field delimiter char (default=\\t)",
                        default = '\t', type = str)
    parser.add_argument("-f", metavar = "FIELD_NUM", dest = "field",
                        help = "target variable field number (starts at 1)",
                        default = 1, type = int)
    parser.add_argument("-n", metavar = "REPEATS", dest = "repeats",
                        help = "statistical test repeats (default=50)",
                        default = 50, type = int)
    # parse CLI
    args = parser.parse_args()
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(file = stderr)
        sys.exit(1)
    input_fn = args.input_fn
    sep_char = args.sep_char
    field = args.field
    repeats = args.repeats
    if field == -1:
        print("-f is mandatory", file = stderr)
        exit(1)
    # compute default params
    all_lines = lines_of_file(input_fn)
    nb_lines = len(all_lines)
    # replace all lines by the values of interest
    all_values = get_all_values(sep_char, field, all_lines)
    # show some stats to the user for troubleshooting
    mini = np.min(all_values)
    aveg = np.average(all_values)
    stdv = np.std(all_values)
    maxi = np.max(all_values)
    print("min,avg+/-std,max: %.3f,%.3f+/-%.3f,%.3f" %
          (mini, aveg, stdv, maxi))
    acc = []
    for i in range(repeats):
        good_vals = []
        dichotomic_search(all_values, good_vals, 0, round(nb_lines / 2), nb_lines)
        best = min(good_vals)
        print("smallest: %d" % best)
        acc.append(best)
    upper_bound = max(acc)
    print("REQUIRED: %d" % upper_bound)

if __name__ == '__main__':
    main()
