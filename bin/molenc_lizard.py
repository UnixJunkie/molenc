#!/usr/bin/env python3

# Copyright (C) 2020, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.
#
# Simple molecular encoder using a handful of molecular descriptors.
# type signature: molecule -> (MolW,cLogP,TPSA,RotB,HBA,HBD,FC)
# Why is this script called molenc_lizard.py
# First, because it is part of the molenc project.
# Second, because if you have a little bit of imagination, a lizard
# is just a very small dragon.
#
# The --remove-aliens option was inspired by
# Tran-Nguyen, V. K., Jacquemard, C., & Rognan, D. (2020).
# "LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening."
# Journal of Chemical Information and Modeling.
# https://doi.org/10.1021/acs.jcim.0c00155

import argparse, rdkit, sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.AtomPairs import Pairs

PeriodicTable = Chem.GetPeriodicTable()

def RobustSmilesMolSupplier(filename):
    with open(filename) as input:
        for i, line in enumerate(input.readlines()):
            words = line.strip().split()
            smile = ""
            name = "NO_NAME"
            try:
                smile = words[0]
                name = words[1]
                yield (i, Chem.MolFromSmiles(smile), name)
            except IndexError:
                # not enough fields on line
                yield (i, None, name)

# detect strange molecules
def is_alien(MolW, cLogP, TPSA, RotB, HBA, HBD, FC):
    # Step 1: TODO?
    #         Organic compound filter. Molecules bearing at
    #         least one atom other than H, C, N, O, P, S, F, Cl, Br, and I
    #         were removed.
    # Step 2: ignored; too complex (using IC50 curve shape + other things).
    # Step 3: Molecular property range filter. Remaining
    #         actives and inactives were kept if:
    # rdkit's MolW unit seems to be g/mol
    # 1C -> 12Da = 12g/mol <=> 1Da = 1g/mol
    return (MolW <= 150 or MolW >= 800 or # 150 < MolW < 800 Da
            cLogP <= -3.0 or cLogP >= 5.0 or # −3.0 < AlogP < 5.0
            RotB >= 15 or # RotB < 15
            HBA >= 10 or # HBA < 10
            HBD >= 10 or # HBD < 10
            FC <= -2 or FC >= 2) # −2.0 < FC < 2.0

# message string describing the alien
def alien_diagnose(i, name, MolW, cLogP, TPSA, RotB, HBA, HBD, FC):
    err_msg = "%d %s" % (i, name)
    if MolW <= 150:
        err_msg += "; (MolW=%g) <= 150" % MolW
    if MolW >= 800:
        err_msg += "; (MolW=%g) >= 800" % MolW
    if cLogP <= -3.0:
        err_msg += "; (cLogP=%g) <= -3" % cLogP
    if cLogP >= 5.0:
        err_msg += "; (cLogP=%g) >= 5" % cLogP
    if RotB >= 15:
        err_msg += "; (RotB=%d) >= 15" % RotB
    if HBA >= 10:
        err_msg += "; (HBA=%d) >= 10" % HBA
    if HBD >= 10:
        err_msg += "; (HBD=%d) >= 10" % HBD
    if FC <= -2:
        err_msg += "; (FC=%d) <= -2" % FC
    if FC >= 2:
        err_msg += "; (FC=%d) >= 2" % FC
    return err_msg

def fun_for_mol_desc(mol_desc):
    if mol_desc == "MolW":
        return Descriptors.MolWt
    elif mol_desc == "HA":
        return Lipinski.HeavyAtomCount
    elif mol_desc == "cLogP":
        return Descriptors.MolLogP
    elif mol_desc == "AR":
        return Lipinski.NumAromaticRings
    elif mol_desc == "TPSA":
        return Descriptors.TPSA
    elif mol_desc == "RotB":
        return Descriptors.NumRotatableBonds
    elif mol_desc == "HBA":
        return Descriptors.NumHAcceptors
    elif mol_desc == "HBD":
        return Descriptors.NumHDonors
    elif mol_desc == "FC":
        return Chem.rdmolops.GetFormalCharge
    else:
        print("FATAL: molenc_lizard.py: fun_for_mol_desc: \
        unsupported: %s" % mol_desc, file=sys.stderr)
        exit(1)

def main():
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Project molecules read from a SMILES file into an 10D \
        space whose dimensions are molecular descriptors: \
        (MolW, HA, cLogP, AR, TPSA, RotB, HBA, HBD, FC)")
    parser.add_argument("-i", metavar = "input_smi", dest = "input_smi",
                        help = "input SMILES file")
    parser.add_argument("-o", metavar = "output_csv", dest = "output_csv",
                        help = "output CSV file")
    parser.add_argument('--no-header', dest='no_header',
                        action='store_true', default=False,
                        help = "no CSV header in output file")
    parser.add_argument('--only', metavar = "mol_desc", dest='mol_desc',
                        default="all",
                        help = "if you want to compute just one molecular \
                        descriptor")
    # just warn about aliens by default
    parser.add_argument('--remove-aliens', dest='rm_aliens',
                        action='store_true', default=False,
                        help = "don't allow aliens in output file")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_smi = args.input_smi
    output_csv = args.output_csv
    rm_aliens = args.rm_aliens
    no_header = args.no_header
    mol_desc = args.mol_desc
    out_count = 0
    alien_count = 0
    error_count = 0
    fun = None
    if mol_desc != "all":
        fun = fun_for_mol_desc(mol_desc)
    with open(output_csv, 'w') as out_file:
        if not no_header:
            if mol_desc == "all":
                print("#name,MolW,HA,cLogP,AR,TPSA,RotB,HBA,HBD,FC",
                      file=out_file)
            else:
                print("#name,%s" % mol_desc, file=out_file)
        for i, mol, name in RobustSmilesMolSupplier(input_smi):
            if mol is None:
                error_count += 1
            else:
                if mol_desc == "all":
                    MolW = Descriptors.MolWt(mol)
                    HA = Lipinski.HeavyAtomCount(mol)
                    cLogP = Descriptors.MolLogP(mol)
                    AR = Lipinski.NumAromaticRings(mol)
                    TPSA = Descriptors.TPSA(mol)
                    RotB = Descriptors.NumRotatableBonds(mol)
                    HBA = Descriptors.NumHAcceptors(mol)
                    HBD = Descriptors.NumHDonors(mol)
                    FC = Chem.rdmolops.GetFormalCharge(mol)
                    alien = is_alien(MolW, cLogP, TPSA, RotB, HBA, HBD, FC)
                    if alien:
                        alien_str = alien_diagnose(i, name, MolW, cLogP, TPSA,
                                                   RotB, HBA, HBD, FC)
                        print("WARN: %s" % alien_str, file=sys.stderr)
                        alien_count += 1
                    if (not alien) or (not rm_aliens):
                        csv_line = "%s,%g,%d,%g,%d,%g,%d,%d,%d,%d" % \
                            (name, MolW, HA, cLogP, AR, TPSA, RotB,
                             HBA, HBD, FC)
                        print(csv_line, file=out_file)
                        out_count += 1
                else:
                    csv_line = "%s,%g" % (name, fun(mol))
                    print(csv_line, file=out_file)
                    out_count += 1
    total_count = out_count + error_count
    if rm_aliens:
        total_count += alien_count
    print("encoded: %d aliens: %d errors: %d total: %d" % \
          (out_count, alien_count, error_count, total_count),
          file=sys.stderr)

if __name__ == '__main__':
    main()
