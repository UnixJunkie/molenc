
import rdkit, deepsmiles, random, re, sys, typing
from rdkit import Chem

def nb_heavy_atom_neighbors(a: rdkit.Chem.rdchem.Atom) -> int:
    res = 0
    for neighb in a.GetNeighbors():
        if neighb.GetAtomicNum() > 1:
            res += 1
    return res

# return (#HA, #H)
def count_neighbors(a: rdkit.Chem.rdchem.Atom) -> tuple[int, int]:
    nb_heavy = nb_heavy_atom_neighbors(a)
    nb_H = a.GetTotalNumHs()
    return (nb_heavy, nb_H)

# DeepSMILES: no rings neither branches opening/closing
to_deep_smiles = deepsmiles.Converter(rings=True, branches=True)

# space-separate all DeepSMILES tokens corresponding to given SMILES
def tokenize_one(smi: str) -> str:
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

def random_reorder_atoms(mol: rdkit.Chem.rdchem.Mol):
    rand_order = list(range(mol.GetNumAtoms()))
    random.shuffle(rand_order)
    rand_mol = Chem.RenumberAtoms(mol, newOrder=rand_order)
    return rand_mol

# return n random versions of smi
def smi_randomize(smi: str, n: int, seed: int) -> list[str]:
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
    def __init__(self, smi: str):
        self.mol = Chem.MolFromSmiles(smi)
        self.mat = Chem.GetDistanceMatrix(self.mol)

    # (atomic_num, #HA, #H, valence - #H, formal_charge)
    def type_atom(self, i: int) -> list[int]:
        a = self.mol.GetAtomWithIdx(i)
        anum = a.GetAtomicNum()
        # do this on the ocaml side, since we get anum
        # assert(anum > 1) # we want to consider only heavy atoms
        nb_HA, nb_H = count_neighbors(a)
        valence = a.GetTotalValence()
        HA_used_val = valence - nb_H
        formal_charge = a.GetFormalCharge()
        return [anum, nb_HA, nb_H, HA_used_val, formal_charge]

    # new (2023) atom typing scheme
    def type_EltFCaroNeighbs(self, i: int) -> list[int]:
        a = self.mol.GetAtomWithIdx(i)
        anum = a.GetAtomicNum()
        fc = a.GetFormalCharge()
        aro = int(a.GetIsAromatic())
        # count direct neighbors
        nb_other = 0 # unsupported atoms
        nb_C  = 0
        nb_H  = a.GetTotalNumHs()
        nb_N  = 0
        nb_O  = 0
        nb_P  = 0
        nb_S  = 0
        nb_F  = 0
        nb_Cl = 0
        nb_Br = 0
        nb_I  = 0
        for neighb in a.GetNeighbors():
            x = neighb.GetAtomicNum()
            if x > 1: # Hs already counted before (including implicits)
                if x == 6:
                    nb_C += 1
                elif x == 7:
                    nb_N += 1
                elif x == 8:
                    nb_O += 1
                elif x == 15:
                    nb_P += 1
                elif x == 16:
                    nb_S += 1
                elif x == 9:
                    nb_F += 1
                elif x == 17:
                    nb_Cl += 1
                elif x == 35:
                    nb_Br += 1
                elif x == 53:
                    nb_I += 1
                else:
                    print("WARN: unsupported anum: %d" % x, file=sys.stderr)
                    nb_other += 1
        return [anum, fc, aro, nb_other, nb_C, nb_H, nb_N, nb_O, nb_P, nb_S, nb_F, nb_Cl, nb_Br, nb_I]

    # simpler atom typing scheme, to reduce the dimension of fingerprints, if needed
    def type_atom_simple(self, i: int) -> list[int]:
        a = self.mol.GetAtomWithIdx(i)
        anum = a.GetAtomicNum()
        fc = a.GetFormalCharge()
        aro = int(a.GetIsAromatic())
        heavies, hydrogens = count_neighbors(a)
        return [anum, fc, aro, heavies, hydrogens]
    
    # # pyml_bindgen doesn't support list of tuples or even tuples...
    # # type each atom of the molecule
    # def type_atoms(self):
    #     res = []
    #     for a in self.mol.GetAtoms():
    #         res.append(type_atom(a))
    #     return res

    def get_num_atoms(self) -> int:
        return self.mol.GetNumAtoms()

    # get the distance (in bonds) between a pair of atoms
    def get_distance(self, i: int, j: int) -> int:
        return int(self.mat[i][j])

    # seed: random_seed
    # n: number of randomized SMILES to use
    # randomize: boolean
    # smi: SMILES to work on
    def get_deep_smiles(self, seed: int, n: int, randomize: bool, smi: str) -> list[str]:
        if n > 1:
            rand_smiles = smi_randomize(smi, n, seed)
            res = list(map(tokenize_one, rand_smiles))
            return res
        else:
            rand_smi = smi
            if randomize:
                rand_smi = smi_randomize(smi, 1, seed)[0]
            return [tokenize_one(rand_smi)]

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
