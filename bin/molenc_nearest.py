#!/usr/bin/env python3

# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, Graduate School of Frontier Sciences,
# The University of Tokyo, Japan.
#
# Annotate a SMILES file using the nearest molecule found in another SMILES
# file according to the Tanimoto score for the given molecular fingerprint
# WARNING: this is O(N^2), so inefficient in case there are many
#          "reference" molecules

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

def parse_smiles_line(line: str) -> tuple[str, str]:
    # print("DEBUG: %s" % line)
    split = line.strip().split()
    smi = split[0]
    name = split[1]
    return (smi, name)

generator = rdFingerprintGenerator.GetMorganGenerator(2, fpSize=2048)

def ecfp4_of_mol(mol):
    return generator.GetFingerprint(mol)

def find_nearest(query_mol, ref_fps):
    best_i = -1
    best_tani = 0.0
    query_fp = ecfp4_of_mol(query_mol)
    for i, ref_fp in enumerate(ref_fps):
        curr_tani = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
        if curr_tani > best_tani:
            best_i = i
            best_tani = curr_tani
    return (best_i, best_tani)

if __name__ == '__main__':
    before = time.time()
    # CLI options parsing
    parser = argparse.ArgumentParser(description = 'find nearest molecule \
    (in fingerprint space)')
    parser.add_argument('-i',
                        metavar = '<filename>', type = str,
                        dest = 'input_fn',
                        help = 'molecules to annotate (SMILES file)')
    parser.add_argument('-r',
                        metavar = '<filename>', type = str,
                        dest = 'ref_fn',
                        default = '',
                        help = 'reference molecules (SMILES file)')
    # parser.add_argument('-np',
    #                     metavar = '<int>', type = int,
    #                     dest = 'nprocs',
    #                     default = 1,
    #                     help = 'max number of processes')
    # parse CLI ---------------------------------------------------------
    if len(sys.argv) == 1:
        # user has no clue of what to do -> usage
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    ref_fn = args.ref_fn
    log("computing FPs for ref. molecules...")
    ref_names = []
    ref_smi = []
    ref_fps = []
    with open(ref_fn, 'r') as input:
        for line in input.readlines():
            strip = line.strip()
            smi, name = parse_smiles_line(strip)
            mol = Chem.MolFromSmiles(smi)
            if mol != None:
                ref_names.append(name)
                ref_smi.append(smi)
                ref_fps.append(ecfp4_of_mol(mol))
    log("iterating over cand. molecules...")
    with open(input_fn, 'r') as input:
        for line in input.readlines():
            strip = line.strip()
            smi, name = parse_smiles_line(strip)
            mol = Chem.MolFromSmiles(smi)
            if mol != None:
                i, best_tani = find_nearest(mol, ref_fps)
                # the candidate molecule
                print('%s\t%s' % (smi, name))
                # its molecular annotation
                print('%s\t%s_T=%.2f' % (ref_smi[i], ref_names[i], best_tani))
    after = time.time()
    dt = after - before
    log('dt: %.2f' % dt)
