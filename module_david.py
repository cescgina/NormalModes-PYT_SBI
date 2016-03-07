### MODULE 1 ####

from Bio.PDB import *
import sys
import re

def isPDBcode(code):
    if len(code)==4:
        return code
    else:
        raise ValueError



def from_pdb_code_to_structure(code):
    """
    From a specific pdb code this function retrieves the structure file
    from the web database and generated a structure instance.
    """
    code=isPDBcode(code)
    pdbl = PDBList()
    parser=PDBParser()
    structure=parser.get_structure(code,pdbl.retrieve_pdb_file(code))
    return structure



def isNMR(structure):
    """
    From a structure instance it returns True if the structure_method is
    'solution nmr'
    """
    if structure.header['structure_method']=='solution nmr':
        return True
    else:
        raise ValueError
