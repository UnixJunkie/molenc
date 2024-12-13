#!/usr/bin/env python3

# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Gaussian Process Regressor CLI wrapper
# input file must have three columns; tab-separated;
# line format: SMILES\tMOLNAME\tpIC50

import argparse
import gpflow
import joblib
import math
import numpy as np
import os
import random
import rdkit
import scipy
import sklearn
import sys
import tempfile
import tensorflow as tf
import time
import typing

from gpflow.mean_functions import Constant
from gpflow.utilities import positive
from gpflow.utilities.ops import broadcasting_elementwise
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy import sparse
from sklearn.metrics import r2_score
# from sklearn.metrics import root_mean_squared_error
from sklearn.metrics import mean_squared_error # if root_mean_squared_error is not available
from sklearn.preprocessing import StandardScaler

# sklearn.metrics.root_mean_squared_error
def root_mean_squared_error(y_true, y_pred):
    return mean_squared_error(y_true, y_pred, squared=False)

# FBR: TODO NxCV should really be parallelized...

# original code from
# https://github.com/Ryan-Rhys/The-Photoswitch-Dataset/blob/master/examples/gp_regression_on_molecules.ipynb
# then refactored by Patrick Walters

class Tanimoto(gpflow.kernels.Kernel):
    def __init__(self):
        super().__init__()
        # constrain kernel variance value to be positive during optimization
        self.variance = gpflow.Parameter(1.0, transform=positive())

    def K(self, X, X2=None):
        """
        Compute the Tanimoto kernel matrix
        σ² * ((<x, y>) / (||x||^2 + ||y||^2 - <x, y>))

        :param X: N x D array
        :param X2: M x D array. If None, compute the N x N kernel matrix for X.
        :return: The kernel matrix of dimension N x M
        """
        if X2 is None:
            X2 = X

        Xs  = tf.reduce_sum(tf.square(X),  axis=-1) # ||X||^2
        X2s = tf.reduce_sum(tf.square(X2), axis=-1) # ||X2||^2
        outer_product = tf.tensordot(X, X2, [[-1], [-1]])

        return self.variance * outer_product / \
            (broadcasting_elementwise(tf.add, Xs, X2s) - outer_product)

    def K_diag(self, X):
        """
        Compute the diagonal of the N x N kernel matrix of X
        :param X: N x D array
        :return: N x 1 array
        """
        return tf.fill(tf.shape(X)[:-1], tf.squeeze(self.variance))

class TanimotoGP:
    def __init__(self, maxiter=100):
        self.m = None
        self.maxiter = maxiter
        self.y_scaler = StandardScaler()

    def objective_closure(self):
        return -self.m.log_marginal_likelihood()

    def fit(self, X_train, y_train):
        y_train_scaled = self.y_scaler.fit_transform(y_train.reshape(-1, 1))
        k = Tanimoto()
        self.m = gpflow.models.GPR(data=(X_train.astype(np.float64), y_train_scaled),
                                   mean_function=Constant(np.mean(y_train_scaled)),
                                   kernel=k,
                                   noise_variance=1)
        opt = gpflow.optimizers.Scipy()
        opt.minimize(self.objective_closure, self.m.trainable_variables,
                     options=dict(maxiter=self.maxiter))

    def predict(self, X_test):
        y_pred, y_var = self.m.predict_f(X_test.astype(np.float64))
        y_pred = self.y_scaler.inverse_transform(y_pred)
        return y_pred.flatten(), y_var.numpy().flatten()

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

def train_test_split(train_portion, lines):
    n = len(lines)
    train_n = math.ceil(train_portion * n)
    test_n = n - train_n
    log('train/test: %d/%d' % (train_n, test_n))
    train = lines[0:train_n]
    test = lines[train_n:]
    assert(len(train) + len(test) == n)
    return (train, test)

def gpr_test(model, X_test):
    preds, _vars = model.predict(X_test)
    return preds

# production predictions; w/ stddev
def gpr_pred(model, X_test):
    # (preds, vars)
    return model.predict(X_test)

def list_take_drop(l, n):
    took = l[0:n]
    dropped = l[n:]
    return (took, dropped)

# cut list into several lists
def list_split(l, n):
    x = float(len(l))
    test_n = math.ceil(x / float(n))
    # log('test_n: %d' % test_n)
    res = []
    taken = []
    rest = l
    for _i in range(0, n):
        curr_test, curr_rest = list_take_drop(rest, test_n)
        curr_train = taken + curr_rest
        res.append((curr_train, curr_test))
        taken = taken + curr_test
        rest = curr_rest
    return res

def parse_smiles_line(line: str) -> tuple[str, str, float]:
    # print("DEBUG: %s" % line)
    split = line.strip().split()
    smi = split[0]
    name = split[1]
    pIC50 = 0.0 # default value
    if len(split) == 3:
        pIC50 = float(split[2])
    return (smi, name, pIC50)

def parse_smiles_lines(lines: list[str]) -> tuple[list[rdkit.Chem.rdchem.Mol],
                                                  list[str],
                                                  list[float]]:
    mols = []
    names = []
    pIC50s = []
    for line in lines:
        smi, name, pIC50 = parse_smiles_line(line)
        mol = Chem.MolFromSmiles(smi)
        if mol != None:
            mols.append(mol)
            names.append(name)
            pIC50s.append(pIC50)
        else:
            log("ERROR: molenc_gpr.py: parse_smiles_lines: could not parse smi for %s: %s" % \
                (name, smi))
    return (mols, names, pIC50s)

def dump_pred_scores(output_fn, names, preds):
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for name, pred in zip(names, preds):
                print('%s\t%f' % (name, pred), file=output)

def dump_pred_vars(output_fn, names, preds, stddevs):
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for name, pred, std in zip(names, preds, stddevs):
                print('%s\t%f\t%f' % (name, pred, std), file=output)

def gnuplot(title0, actual_values, predicted_values):
    # escape underscores so that gnuplot doesn't interprete them
    title = title0.replace('_', '\_')
    xy_min = min(actual_values + predicted_values)
    xy_max = max(actual_values + predicted_values)
    _, commands_temp_fn = tempfile.mkstemp(prefix="gpr_", suffix=".gpl")
    commands_temp_file = open(commands_temp_fn, 'w')
    _, data_temp_fn = tempfile.mkstemp(prefix="gpr_", suffix=".txt")
    data_temp_file = open(data_temp_fn, 'w')
    gnuplot_commands = \
        ["set xlabel 'actual'",
         "set ylabel 'predicted'",
         "set xtics out nomirror",
         "set ytics out nomirror",
         "set xrange [%f:%f]" % (xy_min, xy_max),
         "set yrange [%f:%f]" % (xy_min, xy_max),
         "set key left",
         "set size square",
         "set title '%s'" % title,
         "g(x) = x",
         "f(x) = a*x + b",
         "fit f(x) '%s' u 1:2 via a, b" % data_temp_file.name,
         "plot g(x) t 'perfect' lc rgb 'black', \\",
         "'%s' using 1:2 not, \\" % data_temp_file.name,
         "f(x) t 'fit'"]
    # dump gnuplot commands to temp file
    for l in gnuplot_commands:
        print(l, file=commands_temp_file)
    commands_temp_file.close()
    # dump data to temp file
    for x, y in zip(predicted_values, actual_values):
        print('%f %f' % (x, y), file=data_temp_file)
    data_temp_file.close()
    os.system("gnuplot --persist %s" % commands_temp_fn)

def ecfpX_of_mol(mol: rdkit.Chem.rdchem.Mol, radius) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius, fpSize=2048)
    fp = generator.GetFingerprint(mol)
    arr = np.zeros((1,), int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    # arr: np.ndarray of int64 w/ length 2048
    return arr

def ecfp4_of_mol(mol):
    return ecfpX_of_mol(mol, 2)

# parse SMILES, ignore names, read pIC50s then encode molecules w/ ECFP4 2048b
# return (X_train, y_train)
def read_SMILES_lines_regr(lines):
    mols, names, pIC50s = parse_smiles_lines(lines)
    X_train = np.array([ecfp4_of_mol(mol) for mol in mols])
    y_train = np.array(pIC50s)
    return (X_train, names, y_train)

def gpr_train(X_train, y_train):
    model = TanimotoGP()
    model.fit(X_train, y_train)
    return model

def gpr_train_test_NxCV(all_lines, cv_folds):
    truth = []
    preds = []
    fold = 0
    train_tests = list_split(all_lines, cv_folds)
    for train_set, test_set in train_tests:
        X_train, _names_train, y_train = read_SMILES_lines_regr(train_set)
        X_test, _names_test, y_ref = read_SMILES_lines_regr(test_set)
        model = gpr_train(X_train, y_train)
        truth = truth + list(y_ref)
        y_preds = gpr_test(model, X_test)
        y_preds_lst = list(y_preds)
        r2 = r2_score(y_ref, y_preds_lst)
        rmse = root_mean_squared_error(y_ref, y_preds_lst)
        log('fold: %d R2: %f RMSE: %f' % (fold, r2, rmse))
        preds = preds + y_preds_lst
        fold += 1
    return (truth, preds)

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
    abort_if(model_input_fn != '' and model_output_fn != '',
            "--load and --save are exclusive")
    rng_seed = args.rng_seed
    if rng_seed != -1:
        # only if the user asked for it, we make experiments repeatable
        random.seed(rng_seed)
    output_fn = args.output_fn
    cv_folds = args.cv_folds
    assert(cv_folds >= 1)
    nprocs = args.nprocs
    abort_if(nprocs < 1, 'nprocs must be >= 1')
    train_p = args.train_p
    if model_input_fn != '':
        log('loading trained model: train_p <- 0.0')
        train_p = 0.0
    if model_output_fn != '':
        log('training prod model: train_p <- 1.0')
        train_p = 1.0
    assert(0.0 <= train_p <= 1.0)
    no_compress = args.no_compress
    # work ---------------------------------------------------------
    # read input
    all_lines = lines_of_file(input_fn)
    if train_p > 0.0:
        # if not using a model in production but train-testing:
        # break any ordering the input file might have
        if rng_seed == -1:
            log('WARN: all_lines shuffled w/o rng_seed (not reproducible)')
        random.shuffle(all_lines)
    model = None
    train_lines = []
    test_lines = []
    if cv_folds == 1:
        train_lines, test_lines = train_test_split(train_p, all_lines)
        X_train, _names_train, y_train = read_SMILES_lines_regr(train_lines)
        X_test, names_test, y_test = read_SMILES_lines_regr(test_lines)
        if model_input_fn != '':
            log('loading model from %s' % model_input_fn)
            model = joblib.load(model_input_fn)
        else:
            model = gpr_train(X_train, y_train)
            if model_output_fn != '':
                log('saving model to %s' % model_output_fn)
                if no_compress:
                    joblib.dump(model, model_output_fn, compress=False)
                else:
                    joblib.dump(model, model_output_fn, compress=3)
        # predict w/ trained model
        if train_p == 0.0:
            # production use
            y_preds, y_vars = gpr_pred(model, X_test)
            dump_pred_vars(output_fn, names_test, y_preds, y_vars)
        elif train_p < 1.0:
            y_preds = gpr_test(model, X_test)
            dump_pred_scores(output_fn, names_test, y_preds)
            r2 = r2_score(y_test, y_preds)
            rmse = root_mean_squared_error(y_test, y_preds)
            if train_p > 0.0:
                # train/test case
                log('R2: %.3f RMSE: %.3f fn: %s' % (r2, rmse, input_fn))
            else:
                # maybe production run or predictions
                # on an external validation set
                log('R2: %.3f RMSE: %.3f fn: %s !!! ONLY VALID if test set had target values !!!' % (r2, rmse, input_fn))
    else:
        assert(cv_folds > 1)
        truth, preds = gpr_train_test_NxCV(all_lines, cv_folds)
        log('truths: %d preds: %d' % (len(truth), len(preds)))
        r2 = r2_score(truth, preds)
        rmse = root_mean_squared_error(truth, preds)
        r2_rmse = 'GPR R2=%.3f RMSE=%.3f fn=%s' % (r2, rmse, input_fn)
        log(r2_rmse)
        title = '%s folds=%d %s' % (input_fn, cv_folds, r2_rmse)
        gnuplot(title, truth, preds)
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
