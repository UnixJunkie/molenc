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
import typing
import joblib

from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.metrics import r2_score, root_mean_squared_error
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

# current training-set line format example:
# mol_0,0.0,[1:23;8:3;22:2;45:6;56:1]
def split_AP_line(line):
    fields = line.split(',')
    name = str(fields[0])
    pIC50 = float(fields[1])
    feat_val_strings = fields[2]
    return (name, pIC50, feat_val_strings)

def read_features(row, col, data, i, max_feat_idx, feat_vals_str):
    # print('feat_vals: %d' % (len(feat_vals)))
    feat_vals = feat_vals_str.strip('[]')
    for feat_val in feat_vals.split(';'):
        feat_str, val_str = feat_val.split(':')
        feat = int(feat_str)
        val = float(val_str)
        # print('%d:%d' % (feat, val))
        # assert(feat <= max_feat_idx)        
        if feat > max_feat_idx:
            print("molenc_rfr.py: read_features: \
feat > max_feat_idx: %d > %d" % (feat, max_feat_idx),
                  file=sys.stderr)
            sys.exit(1)
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
        _name, pIC50, features_str = split_AP_line(line)
        # log('%d %f' % (i, pIC50))
        pIC50s.append(pIC50)
        read_features(row, col, data, i, max_feature_index, features_str)
    X = sparse.csr_matrix((data, (row, col)),
                          shape=(nb_lines, max_feature_index + 1))
    log('read_AP_lines_regr: (%d,%d)' % (X.shape[0], X.shape[1]))
    return (X, pIC50s)

def train_test_split(train_portion, lines):
    n = len(lines)
    train_n = math.ceil(train_portion * n)
    test_n = n - train_n
    log('train/test: %d/%d' % (train_n, test_n))
    train = lines[0:train_n]
    test = lines[train_n:]
    assert(len(train) + len(test) == n)
    return (train, test)

def train_ElN_regr(nprocs, X_train, y_train):
    log('ElasticNet regr. training...')
    model = ElasticNetCV(cv=5, n_jobs=nprocs, random_state=314159265,
                         selection='random')
    # model = ElasticNet(alpha=0.004517461806329454, l1_ratio=0.5, random_state=314159265,
    #                    selection='random')
    model.fit(X_train, y_train)
    # e.g. alpha=0.004517461806329454 l1_ratio=0.5
    log('ElasticNet: alpha=%f l1_ratio=%f' % (model.alpha_, model.l1_ratio_))
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

def train_test_NxCV_regr(nprocs, max_features, all_lines, cv_folds):
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
        rmse = root_mean_squared_error(y_ref, y_preds_lst)
        log('fold: %d R2: %f RMSE: %f' % (fold, r2, rmse))
        preds = preds + y_preds_lst
        fold += 1
    return (truth, preds)

if __name__ == '__main__':
    before = time.time()
    train_fn = sys.argv[1]
    test_fn = sys.argv[2]
    max_feat = int(sys.argv[3])
    train_lines = lines_of_file(train_fn)
    test_lines = lines_of_file(test_fn)
    X_train, y_train = read_AP_lines_regr(max_feat, train_lines)
    X_test, y_test = read_AP_lines_regr(max_feat, test_lines)
    model = train_ElN_regr(1, X_train, y_train)
    preds = test(model, X_test)
    log('truths: %d preds: %d' % (len(y_test), len(preds)))
    r2 = r2_score(y_test, preds)
    rmse = root_mean_squared_error(y_test, preds)
    r2_rmse = '%s %s R2=%.3f RMSE=%.4f' % (train_fn, test_fn, r2, rmse)
    log(r2_rmse)
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
