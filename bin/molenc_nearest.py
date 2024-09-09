#!/usr/bin/env python3

# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Annotate a SMILES file using the nearest molecule found in another SMILES
# file according to the Tanimoto score for the given molecular fingerprint
# WARNING: this is O(N^2), so inefficient in case there are many
#          molecules we want to find which one is nearest

import argparse, rdkit, sys, time, typing

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

def hour_min_sec():
    tm = time.localtime()
    return (tm.tm_hour, tm.tm_min, tm.tm_sec)

# use this instead of print to log on stderr
def log(*args, **kwargs):
    hms = hour_min_sec()
    print('%02d:%02d:%02d ' % hms, file=sys.stderr, end='')
    return print(*args, **kwargs, file=sys.stderr)

def abort_if(cond, err_msg):
    if cond:
        log("%s", err_msg)
        sys.exit(1)

def lines_of_file(fn):
    with open(fn) as input:
        return list(map(str.strip, input.readlines()))

def parse_smiles_line(line: str) -> tuple[str, str]:
    # print("DEBUG: %s" % line)
    split = line.strip().split()
    smi = split[0]
    name = split[1]
    return (smi, name)

generator = rdFingerprintGenerator.GetMorganGenerator(2, fpSize=2048)

def ecfp4_of_mol(mol):
    return generator.GetFingerprint(mol)

# parse SMILES, ignore names, read pIC50s then encode molecules w/ ECFP4 2048b
# return (X_train, y_train)
def read_SMILES_lines_regr(lines):
    mols, names, pIC50s = parse_smiles_lines(lines)
    X_train = np.array([ecfp4_of_mol(mol) for mol in mols])
    y_train = np.array(pIC50s)
    return (X_train, names, y_train)

def find_nearest(query_mol, ref_mols):
    best_i = -1
    best_tani = 0.0
    query_fp = ecfp4_of_mol(query_mol)
    for i, mol in enumerate(ref_mols):
        curr_tani = DataStructs.TanimotoSimilarity(query_fp, mol)
        if curr_tani > best_tani:
            best_i = i
            best_tani = curr_tani
    return i

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = 'train/use a GPR model')
    parser.add_argument('-i',
                        metavar = '<filename>', type = str,
                        dest = 'input_fn',
                        help = 'input data file')
    parser.add_argument('-o',
                        metavar = '<filename>', type = str,
                        dest = 'output_fn',
                        default = '',
                        help = 'predictions output file')
    parser.add_argument('--save',
                        metavar = '<filename>', type = str,
                        dest = 'model_output_fn',
                        default = '',
                        help = 'trained model output file')
    parser.add_argument('--load',
                        metavar = '<filename>', type = str,
                        dest = 'model_input_fn',
                        default = '',
                        help = 'trained model input file')
    parser.add_argument('-s',
                        metavar = '<int>', type = int,
                        dest = 'rng_seed',
                        default = -1,
                        help = 'RNG seed')
    parser.add_argument('--no-compress',
                        action = "store_true",
                        dest = 'no_compress',
                        default = False,
                        help = 'turn off saved model compression')
    parser.add_argument('-np',
                        metavar = '<int>', type = int,
                        dest = 'nprocs',
                        default = 1,
                        help = 'max number of processes')
    parser.add_argument('-p',
                        metavar = '<float>', type = float,
                        dest = 'train_p',
                        default = 0.8,
                        help = 'training set proportion')
    parser.add_argument('--NxCV',
                        metavar = '<int>', type = int,
                        dest = 'cv_folds',
                        default = 1,
                        help = 'number of cross validation folds')
    # parse CLI ---------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    model_input_fn = args.model_input_fn
    model_output_fn = args.model_output_fn
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
