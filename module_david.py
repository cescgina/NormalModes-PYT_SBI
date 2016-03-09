"""
This module is a first version of some interestiing function based on biopython
made by DMP.
"""
### MODULE 1 ####

import Bio.PDB as pdb

def pdb_code_check(code):
    """Checks if the code provided have 4 letters"""
    if len(code) == 4:
        return code
    else:
        raise ValueError


def nmr_check(structure):
    """
    From a structure instance it returns True if the structure_method is
    'solution nmr'
    """
    if structure.header['structure_method'] == 'solution nmr':
        return True
    else:
        raise ValueError

def from_pdb_code_to_structure(code):
    """
    From a specific pdb code this function retrieves the structure file
    from the web database and generated a structure instance.
    """
    code = pdb_code_check(code)
    pdbl = pdb.PDBList()
    parser = pdb.PDBParser()
    structure = parser.get_structure(code, pdbl.retrieve_pdb_file(code))
    if nmr_check(structure):
        return structure
    else:
        return ValueError
