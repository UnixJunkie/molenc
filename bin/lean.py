
# automatically find the necessary/minimal training set size using statistics

# find the size of the random sample from the training set that is enough
# to train a model (can help avoiding overfitting, since we might not
# need to use the whole training set to tune the model).
# However, note that for production you will still need to train your
# model on the full training set.

# Bibliography:
# =============
# Domingos, P. (2012). "A few useful things to know about machine learning."
# Communications of the ACM, 55(10), 78-87.

import scipy
from scipy import stats

# initial batch size
b0 = 1000

# subsequent batch increment
bi = 50

# data field separator
sep = ','

# data field
field = 1

# since this is a randomized experiment, we need to run it something like
# 50 times to check the result is somewhat robust
nb_repeats = 50

