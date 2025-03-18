#!/usr/bin/env python3
#
# Extract given tags from an SDF file

import rdkit, sys
from rdkit import Chem

# m.GetProp but w/ a default value
def get_prop_default(m, prop, def_val):
    res = def_val
    try:
        res = m.GetProp(prop)
    except KeyError:
        pass
    return res

def get_props(m, tags):
    res = []
    for t in tags:
        x = get_prop_default(m, t, '')
        res.append(x)
    return res

# FBR: proper CLI
input_fn = sys.argv[1]
output_fn = sys.argv[2]
tags = sys.argv[3].strip().split(',')

with open(output_fn, 'w') as out:
    for mol in Chem.SDMolSupplier(input_fn):
        if mol:
            props = get_props(mol, tags)
            for i, p in enumerate(props):
                if i > 0:
                    print(',%s' % p, end='', file=out)
                else:
                    print('%s' % p, end='', file=out)
            # EOL plus makes empty fields obvious
            print(',', file=out)
