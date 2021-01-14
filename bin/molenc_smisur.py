#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# The "Smiling Surgeon": a surgeon operating directly at the SMILES level

# FRAGMENTATION
#
# 1) cut some cuttable bonds of a molecule
# 2) try to save it as SMILES to see what we get (we should have a mixture)
# 3) atom-type only atoms at the ends of bonds that were cut
# 4) place a dummy atom at each end of the bond that was cut
# 5) each of these dummy atoms should reflect in some way (an isotope number?)
#    the atom type of the previously opposing atom
# 6) save the current state of the edited molecule as a SMILES

# FRAGMENT ASSEMBLY
#
# TODO once and if the former is done
