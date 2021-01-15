
import numpy, rdkit, re

from rdkit import Chem
from enum import IntEnum

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            smile, name = line.strip().split("\t") # enforce TAB-separated
            mol = Chem.MolFromSmiles(smile)
            mol.SetProp("name", name)
            yield (name, mol)

def sort_pairs(pairs):
    res = []
    for (a, b) in pairs:
        x = min(a, b)
        y = max(a, b)
        res.append((x,y))
    return res

# in a bond: atom with lowest index first
# in a list of bonds: bond with lowest first atom index first
def order_bonds_canonically(bonds):
    pairs = map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), bonds)
    min_index_first = sort_pairs(pairs)
    min_index_first.sort()
    return min_index_first

def print_bonds(out, mol):
    print("#bonds:%d" % mol.GetNumBonds(), file=out)
    bonds = order_bonds_canonically(mol.GetBonds())
    for b in bonds:
        print("%d %d" % b, file=out)

def print_distance_matrix(out, mol, threeD):
    if threeD:
        # we use a histogram with bin width 1A
        # this allows to work in 3D, at the cost of much more features
        mat = Chem.Get3DDistanceMatrix(mol)
        diam = int(numpy.max(mat))
        print("#diameter:%d" % diam, file=out)
        nb_atoms = mol.GetNumAtoms()
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                x = int(mat[i][j])
                if j == 0:
                    print("%d" % x, end='', file=out)
                else:
                    print(" %d" % x, end='', file=out)
            print("", file=out) # newline
    else:
        mat = Chem.GetDistanceMatrix(mol)
        diam = numpy.max(mat)
        print("#diameter:%d" % diam, file=out)
        nb_atoms = mol.GetNumAtoms()
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                x = mat[i][j]
                if j == 0:
                    print("%d" % x, end='', file=out)
                else:
                    print(" %d" % x, end='', file=out)
            print("", file=out) # newline

# "#atoms:15 NCGC00261552-01_f00"
def read_atoms_header(line):
    (atoms, nb_atoms, name) = [t(s) for t,s in zip((str,int,str),
                               re.split('[: ]', line))]
    assert(atoms == "#atoms")
    return (nb_atoms, name)

# "0 0,6,2,0"
def read_atom(line):
    (index, nb_pi, atomic_num, nb_HA, charge, stereo) = [t(s) for t,s in
                                                         zip((int,int,int,int,int,int),
                                                         re.split('[, ]', line))]
    return (index, nb_pi, atomic_num, nb_HA, charge, stereo)

# "#bonds:16"
def read_bonds_header(line):
    (bonds, nb_bonds) = [t(s) for t,s in
                         zip((str,int),
                         re.split('[:]', line))]
    assert(bonds == "#bonds")
    return nb_bonds

def bond_type_of_char(c):
    if c == '-':
        return rdkit.Chem.rdchem.BondType.SINGLE
    elif c == ':':
        return rdkit.Chem.rdchem.BondType.AROMATIC
    elif c == '=':
        return rdkit.Chem.rdchem.BondType.DOUBLE
    elif c == '#':
        return rdkit.Chem.rdchem.BondType.TRIPLE
    else:
        assert("molenc_common.py: bond_type_of_char" == "")

def parse_bond_stereo_string(stereo_str):
    if stereo_str == "N":
        return (rdkit.Chem.rdchem.BondStereo.STEREONONE, -1, -1)
    else:
        (stereo_char, a, b) = [t(s) for t,s in zip((str,int,int),
                                                   re.split(':', stereo_str))]
        if stereo_char == "A":
            return (rdkit.Chem.rdchem.BondStereo.STEREOANY, a, b)
        elif stereo_char == "Z":
            return (rdkit.Chem.rdchem.BondStereo.STEREOZ, a, b)
        elif stereo_char == "E":
            return (rdkit.Chem.rdchem.BondStereo.STEREOE, a, b)
        elif stereo_char == "C":
            return (rdkit.Chem.rdchem.BondStereo.STEREOCIS, a, b)
        elif stereo_char == "T":
            return (rdkit.Chem.rdchem.BondStereo.STEREOTRANS, a, b)
        else:
            assert("molenc_common.py: stereo_char not in ['A','Z','E','C','T']" == "")

# ^0 - 1 N$ (no stereo)
# ^1 = 2 Z:0:3$ (some stereo)
def read_bond(line):
    (start_i, c, stop_i, stereo_str) = [t(s) for t,s in zip((int,str,int,str),
                                                            re.split('[ ]', line))]
    stereo = parse_bond_stereo_string(stereo_str)
    return (start_i, bond_type_of_char(c), stop_i, stereo)

class End_of_file(Exception):
    """End of file was reached"""
    pass

# stereo information encoding using integers
class StereoCodes(IntEnum):
    NONE = 0 # default unless specified otherwise
    ANY_CENTER = 1
    R_CENTER = 2
    S_CENTER = 3
    ANY_BOND = 4
    Z_BOND = 5
    E_BOND = 6
    CIS_BOND = 7
    TRANS_BOND = 8

def atom_stereo_code_to_chiral_tag(c):
    if c == 1:
        return Chem.ChiralType.CHI_UNSPECIFIED
    elif c == 2:
        return Chem.ChiralType.CHI_TETRAHEDRAL_CW
    elif c == 3:
        return Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    else:
        assert("molenc_common.py: atom_stereo_code not in [1,2,3]" == "")

def bond_stereo_code_to_bond_stereo(c):
    if c == 0:
        return rdkit.Chem.rdchem.BondStereo.STEREONONE
    elif c == 4:
        return rdkit.Chem.rdchem.BondStereo.STEREOANY
    elif c == 5:
        return rdkit.Chem.rdchem.BondStereo.STEREOZ
    elif c == 6:
        return rdkit.Chem.rdchem.BondStereo.STEREOE
    elif c == 7:
        return rdkit.Chem.rdchem.BondStereo.STEREOCIS
    elif c == 8:
        return rdkit.Chem.rdchem.BondStereo.STEREOTRANS
    else:
        assert("molenc_common.py: bond_stereo_code not in [0,4,5,6,7,8]" == "")

def char_of_bond_stereo(st):
    if st == rdkit.Chem.rdchem.BondStereo.STEREONONE:
        return 'N'
    elif st == rdkit.Chem.rdchem.BondStereo.STEREOANY:
        return 'A'
    elif st == rdkit.Chem.rdchem.BondStereo.STEREOZ:
        return 'Z'
    elif st == rdkit.Chem.rdchem.BondStereo.STEREOE:
        return 'E'
    elif st == rdkit.Chem.rdchem.BondStereo.STEREOCIS:
        return 'C'
    elif st == rdkit.Chem.rdchem.BondStereo.STEREOTRANS:
        return 'T'
    else:
        assert("molenc_common.py: unexpected bond_stereo" == "")
