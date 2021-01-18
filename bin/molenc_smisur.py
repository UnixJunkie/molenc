#!/usr/bin/env python3

# Copyright (C) 2021, Francois Berenger
# Yamanishi laboratory,
# Department of Bioscience and Bioinformatics,
# Faculty of Computer Science and Systems Engineering,
# Kyushu Institute of Technology,
# 680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan.

# The "Smiling Surgeon": a doctor operating directly at the SMILES level

# FRAGMENTATION
#
# 1) OK cut some cuttable bonds of a molecule
#    TODO maybe preserve single bonds coming out of a stereo center
# 2) OK try to save it as SMILES to see what we get (we get a mixture)
# 3) OK atom-type only atoms at the ends of bonds that were cut
# 4) OK use isotope numbers as keys in the (former opposite) atom type map
# 6) OK output this SMILES plus the int->atom_type map as mol_name
#       we should name the fragments also, using parent molecule name + an index
#       # TO RELOAD THIS MAP LATER ON
#       import ast
#       ast.literal_eval(map_str)

# FRAGMENT ASSEMBLY
#
# TODO once and if the former is done
