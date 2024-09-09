#!/usr/bin/env python3

import rdkit, typing
from rdkit import Chem
from rdkit.Chem import AllChem

def create_PEG_chain(length: int):
    s = ''
    for i in range(length):
        s += '[CH2R0X4][CH2R0X4][OH0R0X2]' # SMARTS for one PEG unit
    #return s        
    return Chem.MolFromSmarts(s)

# assume longest linker is 10 units of PEG
peg10 = create_PEG_chain(10)
peg09 = create_PEG_chain(9)
peg08 = create_PEG_chain(8)
peg07 = create_PEG_chain(7)
peg06 = create_PEG_chain(6)
peg05 = create_PEG_chain(5)
peg04 = create_PEG_chain(4)
peg03 = create_PEG_chain(3)
peg02 = create_PEG_chain(2)
#peg01: I assume a single PEG unit is too short to be a proper linker

peg_10_downto_2 = [peg10, peg09, peg08, peg07, peg06, peg05, peg04, peg03, peg02]

# remove the PEG linker, if any
# if not, the molecule is returned unchanged (either it has no linker, or
# the linker is not PEG)
def cut_PEG_linker(mol):
    for patt in peg_10_downto_2:
        if mol.HasSubstructMatch(patt):
            res = AllChem.DeleteSubstructs(mol, patt)
            return (True, res)
    return (False, mol)
