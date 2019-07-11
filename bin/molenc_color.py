#!/usr/bin/env python3

# color molecules from a SMILES file according to per-atom delta score
# values from another file

import matplotlib.pyplot as plot
import rdkit, sys
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor, SimilarityMaps

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split()
            smile = words[0]
            name = " ".join(words[1:]) # everything after the SMILES string
            yield (name, Chem.MolFromSmiles(smile))

# FBR: check params and print usage if necessary

if __name__ == '__main__':
    smiles_fn = sys.argv[1]
    deltas_fn = sys.argv[2]
    delta_max = 0.1 # arbitrary, to normalize deltas and color-scale them
    delta_file = open(deltas_fn, 'r')
    count = 0
    for long_name, mol in RobustSmilesMolSupplier(smiles_fn):
        # split by '_' in case name was postfixed with underscores
        # and additional data
        name = long_name.split('_')[0]
        line = delta_file.readline()
        words = line.split()
        curr_name = words[0]
        if curr_name != name:
            print("names differ: %s != %s" % (name, curr_name))
            exit(1)
        delta_strings = words[1:]
        nb_deltas = len(delta_strings)
        nb_atoms = mol.GetNumAtoms()
        assert(nb_deltas == nb_atoms)
        deltas = list(map(lambda x: float(x), delta_strings))
        rdDepictor.Compute2DCoords(mol) # 2D conformer for figure
        # compute similarity map weights
        weights = []
        for delta in deltas:
            # run-time check that delta is not too high or delta_max too small
            assert(delta <= delta_max)
            weight = delta / delta_max
            weights.append(weight)
        sim_map = Draw.SimilarityMaps.\
                  GetSimilarityMapFromWeights(mol, weights, size = (200,200))
        # the bbox param forces centering the molecule in the figure
        sim_map.savefig(name + '.svg', bbox_inches = 'tight')
        plot.close(sim_map)
        count += 1
    print('processed: %d' % count)
    delta_file.close()
