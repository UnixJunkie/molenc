#!/bin/bash

set -x # DEBUG

diff <(./molenc.py data/caffeine.sdf) <(./molenc.py data/caffeine.smi)
