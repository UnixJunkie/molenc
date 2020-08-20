#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# Short: tool to automatically find the necessary/minimal training set size

# Long: find the size of the random partition from the training set that is enough
# to train a model (can help avoiding overfitting, since we might not
# need to use the whole training set to tune the model).
# However, note that for production you will still need to train your
# model on the full training set.
# In effect, we check that the distribution of the variable to model does not
# significantly change by adding more training samples.
# Usually, with a smaller dataset, you can get a model faster and hence accelerate
# your computational experiments.
#
# Bibliography:
# =============
# Domingos, P. (2012). "A few useful things to know about machine learning."
# Communications of the ACM, 55(10), 78-87.

import argparse, scipy, sys
from scipy import stats

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "automatically determine minimal training set size")
    parser.add_argument("-i", metavar = "INPUT_FILE", dest = "input_fn",
                        help = "line-oriented file containing the \
                        target variable (training-set)")
    parser.add_argument("-d", metavar = "DELIM_CHAR", dest = "sep_char",
                        help = "field delimiter char (default=\\t)", default = '\t')
    parser.add_argument("-f", metavar = "FIELD_NUM", dest = "field",
                        help = "target variable field number (starts from 1)")
    parser.add_argument("-b0", metavar = "INIT_BATCH_SIZE", dest = "init_batch",
                        help = "initial batch size (default=10%% of max)",
                        default = -1)
    parser.add_argument("-b", metavar = "BATCH_SIZE", dest = "batch_incr",
                        help = "batch increment (default=5%% of max)",
                        default = -1)
    # since this is a randomized experiment, we need to run it several
    # times to guarantee some robustness
    parser.add_argument("-n", metavar = "NB_REPEATS", dest = "nb_repeats",
                        help = "number of repetitions (default=50)",
                        default = 50)
    # parse CLI
    args = parser.parse_args()
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(file = sys.stderr)
        sys.exit(1)
    input_fn = args.input_fn
    sep_char = args.sep_char
    field = args.field
    if field == -1:
        print("-f is mandatory", file = sys.stderr)
        exit(1)
    init_batch = args.init_batch
    batch_incr = args.batch_incr
    nb_repeats = args.nb_repeats
    # find default values
    if init_batch == -1:
        exit(1)
    if batch_incr == -1:
        exit(1)
    for i in range(nb_repeats):
        exit(1)

if __name__ == '__main__':
    main()
