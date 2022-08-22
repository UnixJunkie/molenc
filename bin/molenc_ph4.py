#!/usr/bin/env python3
#
# Copyright (C) 2022, Francois Berenger
# Tsuda laboratory, Tokyo University,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# project molecules 3D conformers into the pharmacophore features/points space

import math, os, sys, time
from rdkit import Chem

# ph4 feature SMARTS from the Pharmer software
# definitions from Lidio Meireles and David Ryan Koes
# Article: https://doi.org/10.1021/ci200097m
# Code: https://raw.githubusercontent.com/UnixJunkie/pharmer/master/pharmarec.cpp
aro_smarts = ["a1aaaaa1",
              "a1aaaa1"]

hbd_smarts = ["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
              "[#8!H0&!$([OH][C,S,P]=O)]",
              "[#16!H0]"]

hba_smarts = ["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
	      "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"]

pos_smarts = ["[+,+2,+3,+4]",
              "[$(CC)](=N)N", # amidine
	      "[$(C(N)(N)=N)]", # guanidine
              "[$(n1cc[nH]c1)]"]

neg_smarts = ["[-,-2,-3,-4]",
              "C(=O)[O-,OH,OX1]",
              "[$([S,P](=O)[O-,OH,OX1])]",
	      "c1[nH1]nnn1",
              "c1nn[nH1]n1",
              "C(=O)N[OH1,O-,OX1]",
              "C(=O)N[OH1,O-]",
	      "CO(=N[OH1,O-])",
              "[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]"] # trifluoromethyl sulfonamide

hyd_smarts = ["a1aaaaa1",
	      "a1aaaa1",
	      # branched terminals as one point
	      "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
	      "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      "*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      # simple rings only; need to combine points to get good results for 3d structures
	      "[C&r3]1~[C&r3]~[C&r3]1",
	      "[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
	      "[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
	      "[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
	      "[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
	      "[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
	      # aliphatic chains
	      "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
	      "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
	      "[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
	      # sulfur (apparently)
	      "[$([S]~[#6])&!$(S~[!#6])]"]

def pattern_of_smarts(s):
    return Chem.MolFromSmarts(s)

# compile all SMARTS
aro_patterns = map(pattern_of_smarts, aro_smarts)
hbd_patterns = map(pattern_of_smarts, hbd_smarts)
hba_patterns = map(pattern_of_smarts, hba_smarts)
pos_patterns = map(pattern_of_smarts, pos_smarts)
neg_patterns = map(pattern_of_smarts, neg_smarts)
hyd_patterns = map(pattern_of_smarts, hyd_smarts)

# geometric center of a matched pattern
# WARNING: single-conformer molecule is assumed
def average_match(mol, matched_pattern):
    avg_x = 0.0
    avg_y = 0.0
    avg_z = 0.0
    count = float(len(matched_pattern))
    conf0 = mol.GetConformer(0)
    for i in matched_pattern:
        xyz = conf0.GetAtomPosition(i)
        avg_x += xyz.x
        avg_y += xyz.y
        avg_z += xyz.z
    center = (avg_x / count,
              avg_y / count,
              avg_z / count)
    return center

def find_matches(mol, patterns):
    res = []
    for pat in patterns:
        # get all matches for that pattern
        matched = mol.GetSubstructMatches(pat)
        for m in matched:
            # get the center of each matched group
            avg = average_match(mol, m)
            res.append(avg)
    return res

def find_ARO(mol):
    return find_matches(mol, aro_patterns)

def find_HBD(mol):
    return find_matches(mol, hbd_patterns)

def find_HBA(mol):
    return find_matches(mol, hba_patterns)

def find_POS(mol):
    return find_matches(mol, pos_patterns)

def find_NEG(mol):
    return find_matches(mol, neg_patterns)

def euclid(xyz0, xyz1):
    x0, y0, z0 = xyz0
    x1, y1, z1 = xyz1
    dx = x0 - x1
    dy = y0 - y1
    dz = z0 - z1
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def average(l):
    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0
    n = float(len(l))
    for (x, y, z) in l:
        sum_x += x
        sum_y += y
        sum_z += z
    return (sum_x / n,
            sum_y / n,
            sum_z / n)

cluster_HYD = True # as in Pharmer

def find_HYD(mol):
    hydros = find_matches(mol, hyd_patterns)
    if not cluster_HYD:
        return hydros
    else:
        # regroup all hydrophobic features within 2.0A
        res = []
        n = len(hydros)
        idx2cluster = list(range(n))
        for i in range(n):
            h_i = hydros[i]
            cluster_id = idx2cluster[i]
            for j in range(i+1, n):
                h_j = hydros[j]
                if euclid(h_i, h_j) <= 2.0:
                    # same cluster
                    idx2cluster[j] = cluster_id
        cluster_ids = set(idx2cluster)
        for cid in cluster_ids:
            group = []
            for i, h in enumerate(hydros):
                if idx2cluster[i] == cid:
                    group.append(h)
            res.append(average(group))
        return res

def prfx_print(prfx, x, y, z):
    print("%s %f %f %f" % (prfx, x, y, z))

# FBR: handle CLI options properly

def bild_print(out, color, trans, radius, feats):
    if len(feats) > 0:
        out.write(".color %s\n" % color)
        out.write(".transparency %f\n" % trans)
        for (x, y, z) in feats:
            out.write(".sphere %f %f %f %f\n" % (x, y, z, radius))

def bild_print_ARO(out, feats):
    bild_print(out, "green", 0.75, 1.5, feats)

def bild_print_HYD(out, feats):
    bild_print(out, "grey", 0.75, 1.5, feats)

def bild_print_HBD(out, feats):
    bild_print(out, "white", 0.75, 1.25, feats)

def bild_print_HBA(out, feats):
    bild_print(out, "orange", 0.75, 1.25, feats)

def bild_print_POS(out, feats):
    bild_print(out, "blue", 0.75, 1.0, feats)

def bild_print_NEG(out, feats):
    bild_print(out, "red", 0.75, 1.0, feats)

# better than default readline() never throwing an exception
def read_line_EOF(input):
    line = input.readline()
    if line == "":
        raise EOFError
    else:
        return line

# list all molecule names
def names_of_sdf_file(input_fn):
    res = []
    try:
        with open(input_fn, 'r') as input:
            fst_name = read_line_EOF(input).strip()
            res.append(fst_name)
            while True:
                line = read_line_EOF(input).strip()
                while line != "$$$$":
                    line = read_line_EOF(input).strip()
                next_name = read_line_EOF(input).strip()
                res.append(next_name)
    except EOFError:
        return res

def path_prepend(dir, fn):
    if dir == '':
        return fn
    else:
        return (dir + '/' + fn)

if __name__ == '__main__':
    before = time.time()
    input_fn = sys.argv[1]
    mol_names = names_of_sdf_file(input_fn)
    input_dir = os.path.dirname(input_fn)
    mol_supplier = Chem.SDMolSupplier(input_fn)
    count = 0
    for mol, name in zip(mol_supplier, mol_names):
        aromatics = find_ARO(mol)
        donors = find_HBD(mol)
        acceptors = find_HBA(mol)
        positives = find_POS(mol)
        negatives = find_NEG(mol)
        hydrophobes = find_HYD(mol)
        num_feats = sum(map(len, [aromatics, donors, acceptors, positives, negatives, hydrophobes]))
        print("%d:%s" % (num_feats, name))
        bild_fn = path_prepend(input_dir, name + ".bild")
        bild_out = open(bild_fn, "w")
        for (x, y, z) in aromatics:
            prfx_print("ARO", x, y, z)
        bild_print_ARO(bild_out, aromatics)
        for (x, y, z) in donors:
            prfx_print("HBD", x, y, z)
        bild_print_HBD(bild_out, donors)
        for (x, y, z) in acceptors:
            prfx_print("HBA", x, y, z)
        bild_print_HBA(bild_out, acceptors)
        for (x, y, z) in positives:
            prfx_print("POS", x, y, z)
        bild_print_POS(bild_out, positives)
        for (x, y, z) in negatives:
            prfx_print("NEG", x, y, z)
        bild_print_NEG(bild_out, negatives)
        for (x, y, z) in hydrophobes:
            prfx_print("HYD", x, y, z)
        bild_print_HYD(bild_out, hydrophobes)
        bild_out.close()
        count += 1
    after = time.time()
    dt = after - before
    print("%d molecules at %.2f mol/s" % (count, count / dt), file=sys.stderr)
