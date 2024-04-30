#!/usr/bin/env python3

# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Random Forests Regressor or Classifier (RFR or RFC) CLI wrapper

import argparse
import math
import os
import random
import scipy
import sklearn
import sys
import tempfile
import time
import joblib

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import r2_score, mean_squared_error
from scipy import sparse

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

def count_lines_of_file(fn):
    count = 0
    with open(fn) as input:
        for _l in input:
            count += 1
    return count

def lines_of_file(fn):
    with open(fn) as input:
        return list(map(str.strip, input.readlines()))

# current training-set line format: space-separated
# ^<target_val:float> <feat:int>:<val:float> <feat:int>:<val:float> ...$
def split_AP_line(line):
    fields = line.split(' ')
    pIC50 = float(fields[0])
    feat_val_strings = fields[1:]
    return (pIC50, feat_val_strings)

def read_features(row, col, data, i, max_feat_idx, feat_vals):
    # print('feat_vals: %d' % (len(feat_vals)))
    for feat_val in feat_vals:
        feat_str, val_str = feat_val.split(':')
        feat = int(feat_str)
        val = float(val_str)
        # print('%d:%d' % (feat, val))
        assert(feat <= max_feat_idx)
        if val != 0.0:
            row.append(i)
            col.append(feat)
            data.append(val)

def read_AP_lines_regr(max_feature_index, lines):
    nb_lines = len(lines)
    pIC50s = []
    row = []
    col = []
    data = []
    for i, line in enumerate(lines):
        pIC50, features_str = split_AP_line(line)
        # log('%d %f' % (i, pIC50))
        pIC50s.append(pIC50)
        read_features(row, col, data, i, max_feature_index, features_str)
    X = sparse.csr_matrix((data, (row, col)),
                          shape=(nb_lines, max_feature_index + 1))
    log('read_AP_lines_regr: (%d,%d)' % (X.shape[0], X.shape[1]))
    return (X, pIC50s)

def label_of_name(name):
    if name.startswith('active'):
        return 1
    else:
        return 0

def read_AP_lines_class(max_feature_index, lines):
    nb_lines = len(lines)
    names = []
    labels = []
    row = []
    col = []
    data = []
    for i, line in enumerate(lines):
        name, _pIC50, features_str = split_AP_line(line)
        names.append(name)
        labels.append(label_of_name(name))
        read_features(row, col, data, i, max_feature_index, features_str)
    X = sparse.csr_matrix((data, (row, col)),
                          shape=(nb_lines, max_feature_index + 1))
    log('read_AP_lines_class: (%d,%d)' % (X.shape[0], X.shape[1]))
    return (names, X, labels)

def train_test_split(train_portion, lines):
    n = len(lines)
    train_n = math.ceil(train_portion * n)
    test_n = n - train_n
    log('train/test: %d/%d' % (train_n, test_n))
    train = lines[0:train_n]
    test = lines[train_n:]
    assert(len(train) + len(test) == n)
    return (train, test)

def train_regr(ntrees, crit, nprocs, msl, max_features, max_samples,
               X_train, y_train):
    log('regressor training...')
    model = RandomForestRegressor(n_estimators = ntrees,
                                  criterion = crit,
                                  n_jobs = nprocs,
                                  oob_score = True,
                                  min_samples_leaf = msl,
                                  max_features = max_features,
                                  max_samples = max_samples)
    model.fit(X_train, y_train)
    return model

# WARNING: in binary classification, y_train must be only 0s and 1s
def train_class(ntrees, nprocs, msl, max_features, X_train, y_train):
    log('classifier training...')
    model = RandomForestClassifier(n_estimators = ntrees,
                                   criterion = 'gini',
                                   n_jobs = nprocs,
                                   oob_score = True,
                                   min_samples_leaf = msl,
                                   max_features = max_features)
    model.fit(X_train, y_train)
    return model

def test(model, X_test):
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

def train_test_NxCV_regr(max_feat, ntrees, crit, nprocs, msl, max_features,
                         max_samples, all_lines, cv_folds):
    truth = []
    preds = []
    fold = 0
    train_tests = list_split(all_lines, cv_folds)
    for train_set, test_set in train_tests:
        X_train, y_train = read_AP_lines_regr(max_feat, train_set)
        model = train_regr(ntrees, crit, nprocs, msl, max_features,
                           max_samples, X_train, y_train)
        X_test, y_ref = read_AP_lines_regr(max_feat, test_set)
        truth = truth + y_ref
        y_preds = test(model, X_test)
        y_preds_lst = []
        for y in y_preds:
            y_preds_lst.append(y)
        r2 = r2_score(y_ref, y_preds_lst)
        rmse = mean_squared_error(y_ref, y_preds_lst, squared=False)
        log('fold: %d R2: %f RMSE: %f' % (fold, r2, rmse))
        preds = preds + y_preds_lst
        fold += 1
    return (truth, preds)

def train_test_NxCV_class(max_feat, ntrees, nprocs, msl, max_features,
                          all_lines, cv_folds):
    truth = []
    preds = []
    fold = 0
    train_tests = list_split(all_lines, cv_folds)
    for train_set, test_set in train_tests:
        _names, X_train, y_train = read_AP_lines_class(max_feat, train_set)
        model = train_class(ntrees, nprocs, msl, max_features,
                            X_train, y_train)
        _names, X_test, y_ref = read_AP_lines_class(max_feat, test_set)
        truth = truth + y_ref
        y_preds = test(model, X_test)
        y_preds_lst = []
        for y in y_preds:
            y_preds_lst.append(y)
        auc = sklearn.metrics.roc_auc_score(y_ref, y_preds_lst)
        mcc = sklearn.metrics.matthews_corrcoef(y_ref, y_preds_lst)
        log('fold: %d AUC: %f MCC: %.3f' % (fold, auc, mcc))
        preds = preds + y_preds_lst
        fold += 1
    return (truth, preds)

def dump_pred_scores(output_fn, preds):
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for pred in preds:
                print('%f' % pred, file=output)

# a prediction is an integer class label, not a float
# FBR: maybe output the active class proba as an additional column
def dump_pred_labels(output_fn, names, preds):
    assert(len(names) == len(preds))
    if output_fn != '':
        with open(output_fn, 'w') as output:
            for name_pred in zip(names, preds):
                print('%s\t%d' % name_pred, file=output)

def gnuplot(title0, actual_values, predicted_values):
    # escape underscores so that gnuplot doesn't interprete them
    title = title0.replace('_', '\_')
    xy_min = min(actual_values + predicted_values)
    xy_max = max(actual_values + predicted_values)
    _, commands_temp_fn = tempfile.mkstemp(prefix="pcp_rfr_", suffix=".gpl")
    commands_temp_file = open(commands_temp_fn, 'w')
    _, data_temp_fn = tempfile.mkstemp(prefix="pcp_rfr_", suffix=".txt")
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

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = 'train/use a RFR model')
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
    parser.add_argument('-n',
                        metavar = '<int>', type = int,
                        dest = 'ntrees',
                        default = 100,
                        help = 'number of trees')
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
    parser.add_argument('-msl',
                        metavar = '<int>', type = int,
                        dest = 'msl',
                        default = 1,
                        help = 'min samples leaf')
    parser.add_argument('-p',
                        metavar = '<float>', type = float,
                        dest = 'train_p',
                        default = 0.8,
                        help = 'training set proportion')
    parser.add_argument('-c',
                        metavar = '<string>', type = str,
                        dest = 'criterion',
                        default = 'squared_error',
                        help = 'mae|squared_error')
    parser.add_argument('--mtry',
                        metavar = '<string>', type = float,
                        dest = 'train_feats',
                        default = 1.0,
                        help = 'float in ]0,1.0]')
    parser.add_argument('-m',
                        metavar = '<int>', type = int,
                        dest = 'max_feat',
                        default = 10_000,
                        help = 'max feature index (cf. AP dictionary)')
    parser.add_argument('--NxCV',
                        metavar = '<int>', type = int,
                        dest = 'cv_folds',
                        default = 1,
                        help = 'number of cross validation folds')
    parser.add_argument('-ms',
                        metavar = '<float>', type = float,
                        dest = 'max_samples',
                        # keep None here, not 1.0; reported sklearn bug
                        default = None,
                        help = 'max samples per bootstrap')
    parser.add_argument('--classif', dest = 'classification',
                        action ='store_true',
                        default = False,
                        help = 'train a classifier \
                        (default: train a regressor)')
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
    max_feat = args.max_feat # feat. dict related; _NOT_ model hyper param
    # model's hyper parameters
    cv_folds = args.cv_folds
    assert(cv_folds >= 1)
    msl = args.msl
    crit = args.criterion
    assert(crit in ['mae','squared_error'])
    # REMARK: faster to train w/ squared_error, but better model with mae
    ntrees = args.ntrees
    abort_if(ntrees < 50, 'ntrees must be >= 50')
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
    max_features = args.train_feats
    assert(0.0 < max_features <= 1.0)
    # FBR: for classification and regression, the max_features defaults should not be the same
    #      IIRC for classif. this is sqrt(N); for regre. this is N
    max_samples = args.max_samples # this one must default to None
    if max_samples == 1.0:
        max_samples = None # BUG in sklearn RFR probably; this forces 1.0
    no_compress = args.no_compress
    classification = args.classification
    regressor = not classification
    # work ---------------------------------------------------------
    # read in a .AP file w/ pIC50 values
    all_lines = lines_of_file(input_fn)
    if train_p > 0.0:
        # if not using a model in production but train-testing:
        # break any ordering the input file might have
        random.shuffle(all_lines)
    model = None
    train_lines = []
    test_lines = []
    if cv_folds == 1:
        train_lines, test_lines = train_test_split(train_p, all_lines)
        if classification:
            _names_train, X_train, y_train = \
                read_AP_lines_class(max_feat, train_lines)
            names_test, X_test, y_test = \
                read_AP_lines_class(max_feat, test_lines)
            if model_input_fn != '':
                log('load model from %s' % model_input_fn)
                model = joblib.load(model_input_fn)
            else:
                model = train_class(ntrees, nprocs, msl, max_features,
                                    X_train, y_train)
                if model_output_fn != '':
                    log('saving model to %s' % model_output_fn)
                    joblib.dump(model, model_output_fn, compress=('xz',9))
            # predict w/ trained model
            if train_p < 1.0:
                y_preds = test(model, X_test)
                dump_pred_labels(output_fn, names_test, y_preds)
            # train_p = 0.0: we are in production,
            # assume there are no labels
            if train_p > 0.0:
                auc = sklearn.metrics.roc_auc_score(y_test, y_preds)
                mcc = sklearn.metrics.matthews_corrcoef(y_test, y_preds)
                log('AUC: %.3f MCC: %.3f' % (auc, mcc))
        else: #regression
                X_train, y_train = read_AP_lines_regr(max_feat, train_lines)
                X_test, y_test = read_AP_lines_regr(max_feat, test_lines)
                if model_input_fn != '':
                    log('loading model from %s' % model_input_fn)
                    model = joblib.load(model_input_fn)
                else:
                    model = train_regr(ntrees, crit, nprocs, msl, max_features,
                                       max_samples, X_train, y_train)
                    if model_output_fn != '':
                        log('saving model to %s' % model_output_fn)
                        if no_compress:
                            joblib.dump(model, model_output_fn, compress=False)
                        else:
                            joblib.dump(model, model_output_fn, compress=3)
                # predict w/ trained model
                if train_p < 1.0:
                    y_preds = test(model, X_test)
                    dump_pred_scores(output_fn, y_preds)
                    r2 = r2_score(y_test, y_preds)
                    rmse = mean_squared_error(y_test, y_preds, squared=False)
                    if train_p > 0.0:
                        # train/test case
                        log('R2: %f RMSE: %f' % (r2, rmse))
                    else:
                        # maybe production run or predictions
                        # on an external validation set
                        log('R2: %f RMSE: %f !!! ONLY VALID if test set had target values !!!' % (r2, rmse))
    else:
        assert(cv_folds > 1)
        if classification:
                truth, preds = train_test_NxCV_class(max_feat, ntrees, nprocs, msl,
                                                     max_features, all_lines,
                                                     cv_folds)
                log('truths: %d preds: %d' % (len(truth), len(preds)))
                auc = sklearn.metrics.roc_auc_score(truth, preds)
                mcc = sklearn.metrics.matthews_corrcoef(truth, preds)
                log('AUC: %.3f MCC: %.3f' % (auc, mcc))
        else: #regression
                truth, preds = train_test_NxCV_regr(max_feat, ntrees, crit,
                                                    nprocs, msl, max_features,
                                                    max_samples, all_lines,
                                                    cv_folds)
                log('truths: %d preds: %d' % (len(truth), len(preds)))
                r2 = r2_score(truth, preds)
                rmse = mean_squared_error(truth, preds, squared=False)
                r2_rmse = 'R2=%.3f RMSE=%.4f' % (r2, rmse)
                log(r2_rmse)
                title = '%s N=%d folds=%d mtry=%g %s' % (input_fn, ntrees, cv_folds, max_features, r2_rmse)
                gnuplot(title, truth, preds)
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
