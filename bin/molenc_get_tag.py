#!/usr/bin/env python3
#
# Extract given tags from an SDF file
# could also be called sdf2csv

import argparse, rdkit, sys
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

if __name__ == '__main__':
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "compute atom types and distances")
    parser.add_argument("-i", metavar = "input.sdf", dest = "input_fn",
                        help = "molecules input file")
    parser.add_argument("-o", metavar = "output.txt", dest = "output_fn",
                        help = "output file")
    parser.add_argument('-t', metavar = 'tag1,tag2,...', dest = "tags",
                        help = "comma-separated list of tags to extract")
    # parse CLI
    if len(sys.argv) == 1:
        # show help in case user has no clue of what to do
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_fn = args.input_fn
    output_fn = args.output_fn
    tags = args.tags
    # -------------------------------------------------------------------------
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
