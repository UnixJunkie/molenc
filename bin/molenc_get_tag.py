#!/usr/bin/env python3
#
# Extract given tags from an SDF file; not parsing molecules for speed

import argparse, gzip, re, sys

class End_of_file(Exception):
    pass

def open_maybe_compressed_file(fn):
    if fn.endswith(".gz"):
        return gzip.open(fn, 'rt')
    else:
        return open(fn, 'rt')

# regexp to capture an SDF tag (key)
# next line must contain the corresponding value
sdf_tag_name = re.compile('>  <(.+)>')

if __name__ == '__main__':
    # CLI options parsing
    parser = argparse.ArgumentParser(
        description = "Extract listed properties from SDF file")
    parser.add_argument("-i", metavar = "input.sdf[.gz]", dest = "input_fn",
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
    tags = args.tags.strip().split(',')
    tags_set = set(tags)
    # -------------------------------------------------------------------------
    res: dict[str,list[str]] = {}
    for tag in tags:
        res[tag] = []
    # capture all tags we are interested in
    with open_maybe_compressed_file(input_fn) as input:
        try:
            while True:
                line = input.readline()
                if line == '':
                    raise End_of_file
                elif line[0] == '>':
                    match = sdf_tag_name.match(line)
                    if match:
                        key = match.group(1)
                        if key in tags_set:
                            value: str = input.readline().strip()
                            prev_vals = res[key]
                            prev_vals.append(value)
                            res[key] = prev_vals
        except End_of_file:
            pass
    # check all tags of interest were seen the same number of times
    expected = len(res[tags[0]])
    for tag in tags:
        if expected != len(res[tag]):
            print('FATAL: expected: %d <> len(res[%s]) = %d' % \
                  (expected, tag, len(res[tag])), file=sys.stderr)
            exit(1)
    # CSV output
    values = []
    for tag in tags:
        values.append(res[tag])
    with open(output_fn, 'w') as csv_out:
        for i in range(expected):
            for j in range(len(values)):
                print("%s" % values[j][i], file=csv_out, end=',')
            print("", file=csv_out) # EOL
