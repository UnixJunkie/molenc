
import rdkit, deepsmiles, random, re
from rdkit import Chem

def nb_heavy_atom_neighbors(a):
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() > 1:
            res += 1
    return res

# return (#HA, #H)
def count_neighbors(a):
    nb_heavy = nb_heavy_atom_neighbors(a)
    nb_H = a.GetTotalNumHs()
    return (nb_heavy, nb_H)

# DeepSMILES: no rings neither branches opening/closing
to_deep_smiles = deepsmiles.Converter(rings=True, branches=True)

# space-separate all DeepSMILES tokens corresponding to given SMILES
def tokenize_one(smi):
    assert(smi.find('.') == -1) # enforce standardization/salt removal
    mol = Chem.MolFromSmiles(smi)
    # don't canonicalize: the input SMILES might have been randomized on purpose
    protected_smi = Chem.MolToSmiles(mol, allHsExplicit=True, canonical=False)
    protected_dsmi = to_deep_smiles.encode(protected_smi)
    # print("pdsmi: '%s'" % protected_dsmi)
    # space before [ and after ]
    pdsmi = re.sub(r"(\[[^\]]+\])", r" \1 ", protected_dsmi)
    # space before %
    pdsmi = pdsmi.replace('%', ' %')
    # protect branch closings (no branch openings in DeepSMILES)
    pdsmi = pdsmi.replace(')', ' ) ')
    # protect bonds
    pdsmi = pdsmi.replace('-', ' - ')
    # protect - when it is a formal charge
    pdsmi = re.sub(' - (\d)\]', '-\1]', pdsmi)
    pdsmi = re.sub(' - \]', '-]', pdsmi)
    pdsmi = pdsmi.replace('=', ' = ')
    pdsmi = pdsmi.replace('#', ' # ')
    pdsmi = pdsmi.replace('$', ' $ ')
    pdsmi = pdsmi.replace(':', ' : ')
    # protect long numbers (prefixed by % in SMILES)
    pdsmi = re.sub(r"%(\d)(\d)", r" %\1\2 ", pdsmi)
    # single digit numbers are separate words
    pdsmi = re.sub(r" (\d)(\d)", r" \1 \2", pdsmi)
    # protect stereo bonds
    pdsmi = pdsmi.replace("/", " / ")
    pdsmi = pdsmi.replace("\\", " \\ ")
    # several spaces to one
    pdsmi = re.sub('[ ]+', ' ', pdsmi)
    # rm leading/trailing whitespaces
    pdsmi = pdsmi.strip()
    # print("pdsmi: '%s'" % pdsmi)
    return pdsmi

def random_reorder_atoms(mol):
    rand_order = list(range(mol.GetNumAtoms()))
    random.shuffle(rand_order)
    rand_mol = Chem.RenumberAtoms(mol, newOrder=rand_order)
    return rand_mol

# return n random versions of smi
def smi_randomize(smi, n, seed):
    res = []
    mol = Chem.MolFromSmiles(smi)
    random.seed(seed)
    for i in range(n):
        rand_mol = random_reorder_atoms(mol)
        rand_smi = Chem.MolToSmiles(rand_mol, canonical=False)
        res.append(rand_smi)
    return res

class Rdkit:
    # this is needed because the OCaml side want to know how
    # to get an object of type t
    def __init__(self, smi):
        self.mol = Chem.MolFromSmiles(smi)
        self.mat = Chem.GetDistanceMatrix(self.mol)

    # (atomic_num, #HA, #H, valence - #H, formal_charge)
    def type_atom(self, i):
        a = self.mol.GetAtomWithIdx(i)
        anum = a.GetAtomicNum()
        # do this on the ocaml side, since we get anum
        # assert(anum > 1) # we want to consider only heavy atoms
        nb_HA, nb_H = count_neighbors(a)
        valence = a.GetTotalValence()
        HA_used_val = valence - nb_H
        formal_charge = a.GetFormalCharge()
        return [anum, nb_HA, nb_H, HA_used_val, formal_charge]

    # # pyml_bindgen doesn't support list of tuples or even tuples...
    # # type each atom of the molecule
    # def type_atoms(self):
    #     res = []
    #     for a in self.mol.GetAtoms():
    #         res.append(type_atom(a))
    #     return res

    def get_num_atoms(self):
        return self.mol.GetNumAtoms()

    # get the distance (in bonds) between a pair of atoms
    def get_distance(self, i, j):
        return int(self.mat[i][j])

# # tests
# m = Chem.MolFromSmiles('c1ccccc1')
# assert(get_distances(m) == [(0, 1, 1),
#                             (0, 2, 2),
#                             (0, 3, 3),
#                             (0, 4, 2),
#                             (0, 5, 1),
#                             (1, 2, 1),
#                             (1, 3, 2),
#                             (1, 4, 3),
#                             (1, 5, 2),
#                             (2, 3, 1),
#                             (2, 4, 2),
#                             (2, 5, 3),
#                             (3, 4, 1),
#                             (3, 5, 2),
#                             (4, 5, 1)])
# assert(type_atoms(m) == [(6, 2, 1, 3, 0),
#                          (6, 2, 1, 3, 0),
#                          (6, 2, 1, 3, 0),
#                          (6, 2, 1, 3, 0),
#                          (6, 2, 1, 3, 0),
#                          (6, 2, 1, 3, 0)])
