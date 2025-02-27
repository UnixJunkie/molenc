#!/usr/bin/env python3

# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Gaussian Process Classifier CLI wrapper
# input line format: '^SMILES\t[active]MOLNAME$'

import argparse
import joblib
import math
import numpy as np
import os
import random
import rdkit
import sklearn
import sys
import tempfile
import time
import typing

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.metrics import roc_auc_score

def hour_min_sec() -> tuple[float, float, float]:
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

def predict_classes(model, X_test) -> np.ndarray:
    return model.predict(X_test)

def predict_probas(model, X_test) -> np.ndarray:
    probas = model.predict_proba(X_test)
    # let's hope probas[:,1] are probabilities for the actives' class
    return probas[:,1]

def list_take_drop(l: list, n: int) -> tuple[list, list]:
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

def parse_smiles_line(line: str) -> tuple[str, str, bool]:
    # print("DEBUG: %s" % line)
    split = line.strip().split()
    assert(len(split) == 2) # ^SMILES\t[active]name$'
    smi = split[0]
    name = split[1]
    label = name.startswith('active')
    return (smi, name, label)

def parse_smiles_lines(lines: list[str]) -> tuple[list[rdkit.Chem.rdchem.Mol],
                                                  list[str],
                                                  list[bool]]:
    mols = []
    names = []
    labels = []
    for line in lines:
        smi, name, label = parse_smiles_line(line)
        mol = Chem.MolFromSmiles(smi)
        if mol != None:
            mols.append(mol)
            names.append(name)
            labels.append(label)
        else:
            log("ERROR: molenc_gpc.py: parse_smiles_lines: could not parse smi for %s: %s" % \
                (name, smi))
    return (mols, names, labels)

def dump_pred_probas(output_fn, names, probas):
    if len(names) != len(probas):
        log('|names|=%d <> |probas|=%d' % (len(names), len(probas)))
        exit(1)
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for name, p in zip(names, probas):
                print('%s\t%f' % (name, p), file=output)

def dump_pred_labels(output_fn, names, labels):
    if len(names) != len(labels):
        log('|names|=%d <> |probas|=%d' % (len(names), len(labels)))
        exit(1)
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for name, label in zip(names, labels):
                print('%s\t%d' % (name, label), file=output)

def dump_score_labels(output_fn: str,
                      probas: list[float],
                      labels: list[bool]):
    num_p = len(probas)
    num_l = len(labels)
    if num_p != num_l:
        log('|probas|=%d <> |labels|=%d' % (num_p, num_l))
        exit(1)
    with open(output_fn, 'w') as output:
        for score, label in zip(probas, labels):
            print('%f\t%d' % (score, label), file=output)

# use the croc-curve command to compute the ROC curve then return its AUC
def roc_curve(in_score_labels_fn: str,
              out_curve_fn: str) -> float:
    out_auc_fn = out_curve_fn + ".auc"
    cmd = "croc-curve < %s >%s 2>%s" % (in_score_labels_fn, out_curve_fn, out_auc_fn)
    ret = os.system(cmd)
    assert(ret == 0)
    auc = float('nan')
    with open(out_auc_fn) as input:
        line = input.readline()
        strip = line.strip()
        tokens = strip.split()
        auc = float(tokens[4])
    return auc

def temp_file(prfx, sfx):
    _, temp_fn = tempfile.mkstemp(prefix=prfx, suffix=sfx)
    return temp_fn

def gnuplot(title0, auc_curve_fn):
    # escape underscores so that gnuplot doesn't interprete them
    title = title0.replace('_', '\_')
    gnuplot_commands = \
        ["set xlabel 'FPR'",
         "set ylabel 'TPR'",
         "set tics out nomirror",
         "set size square",
         "set xrange [0:1]",
         "set yrange [0:1]",
         "set key left",
         "diag(x) = x",
         "set title '%s'" % title,
         "plot diag(x) not lc rgb 'black', '%s' u 1:2 w l t 'ROC'" % auc_curve_fn]
    # dump gnuplot commands
    commands_temp_fn = temp_file("gpr_", ".gpl")    
    with open(commands_temp_fn, 'w') as output:
        for l in gnuplot_commands:
            print(l, file=output)
    os.system("gnuplot --persist %s" % commands_temp_fn)
    os.remove(commands_temp_fn) # cleanup

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
def read_SMILES_lines_class(lines):
    mols, names, labels = parse_smiles_lines(lines)
    X_train = np.array([ecfp4_of_mol(mol) for mol in mols])
    y_train = np.array(labels)
    return (X_train, names, y_train)

def gpc_train(X_train, y_train, seed=0):
    if len(X_train) != len(y_train):
        log('|X_train|=%d <> |y_train|=%d' % (len(X_train), len(y_train)))
        exit(1)
    rbf_k = 1.0 * RBF(1.0)
    gpc = GaussianProcessClassifier(kernel=rbf_k,
                                    random_state=seed,
                                    optimizer='fmin_l_bfgs_b',
                                    n_restarts_optimizer=5,
                                    copy_X_train=True)
    gpc.fit(X_train, y_train)
    return gpc

def gpc_train_test_NxCV(all_lines, cv_folds):
    truth = []
    proba_preds = []
    fold = 0
    train_tests = list_split(all_lines, cv_folds)
    for train_set, test_set in train_tests:
        X_train, _names_train, y_train = read_SMILES_lines_class(train_set)
        X_test, _names_test, y_ref = read_SMILES_lines_class(test_set)
        model = gpc_train(X_train, y_train)
        truth = truth + list(y_ref)
        pred_probas = predict_probas(model, X_test)
        roc_auc = roc_auc_score(y_ref, pred_probas)
        log('fold: %d AUC: %.3f' % (fold, roc_auc))
        proba_preds = proba_preds + list(pred_probas)
        fold += 1
    return (truth, proba_preds)

def show_roc_curve(plot_title, pred_probas, true_labels):
    score_labels_fn = temp_file("gpc_", ".score_labels")
    auc_curve_fn = temp_file("gpc_", ".roc")
    dump_score_labels(score_labels_fn, pred_probas, true_labels)
    _croc_curve_auc = roc_curve(score_labels_fn, auc_curve_fn)
    # print('DEBUG: ROC AUC file: %s' % auc_curve_fn)
    gnuplot(plot_title, auc_curve_fn)
    # cleanup
    os.remove(score_labels_fn)
    os.remove(auc_curve_fn)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = 'train/use a GPC model')
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
    parser.add_argument('--no-plot',
                        action = "store_true",
                        dest = 'no_plot',
                        default = False,
                        help = 'do not show the ROC curve')
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
        # user has no clue of what to do: usage
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
    no_plot = args.no_plot
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
        X_train, names_train, y_train = read_SMILES_lines_class(train_lines)
        X_test, names_test, y_test = read_SMILES_lines_class(test_lines)
        assert(len(X_train) == len(names_train) == len(y_train))
        assert(len(X_test) == len(names_test) == len(y_test))
        if model_input_fn != '':
            log('loading model from %s' % model_input_fn)
            model = joblib.load(model_input_fn)
        else:
            # print('|X_train|=%d' % len(X_train))
            # print('|y_train|=%d' % len(y_train))
            model = gpc_train(X_train, y_train)
            if model_output_fn != '':
                log('saving model to %s' % model_output_fn)
                if no_compress:
                    joblib.dump(model, model_output_fn, compress=False)
                else:
                    joblib.dump(model, model_output_fn, compress=3)
        # predict w/ trained model
        if train_p < 1.0:
            pred_probas = predict_probas(model, X_test)
            # print('|X_test|=%d' % len(X_test))
            # print('|y_test|=%d' % len(y_test))
            # print('|pred_probas|=%d' % len(pred_probas))
            dump_pred_probas(output_fn, names_test, pred_probas)
            # print('t(y_test)=%s' % type(y_test))
            # print('t(pred_probas)=%s' % type(pred_probas))
            auc = roc_auc_score(y_test, pred_probas)
            if train_p > 0.0:
                # train/test case
                title = "GPC AUC=%.3f fn=%s" % (auc, input_fn)
                log(title)
                if not no_plot:
                    show_roc_curve(title, pred_probas, y_test)
            else:
                # maybe production run or predictions
                # on an external validation set
                log('AUC: %.3f fn: %s !!! ONLY VALID if test set had target values !!!' %
                    (auc, input_fn))
    else:
        assert(cv_folds > 1)
        true_labels, preds = gpc_train_test_NxCV(all_lines, cv_folds)
        assert(len(true_labels) == len(preds))
        auc = roc_auc_score(true_labels, preds)
        auc_msg = 'GPC folds=%d AUC=%.3f fn=%s' % (cv_folds, auc, input_fn)
        log(auc_msg)
        # show the overall (all folds combined) ROC AUC curve
        if not no_plot:
            show_roc_curve(auc_msg, preds, true_labels)
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
