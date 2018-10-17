
# in a bond: atom with lowest index first
# in a list of bonds: bond with lowest first atom index first
def order_bonds_canonically(bonds):
    pairs = map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), bonds)
    min_index_first = map(lambda (a, b): (min(a, b), max(a, b)), pairs)
    min_index_first.sort()
    return min_index_first
