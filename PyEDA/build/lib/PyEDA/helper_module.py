"""
This module is a first version of some interestiing function based on biopython
made by DMP.
"""
### MODULE 1 ####

import urllib
import sys
import Bio.PDB as pdb



def store_header_text(filename):
    """It gets the header of the structure from the ent file."""
    try:
        file_hand = open(filename, "r")
        current_header = ''
        for line in file_hand:
            if line.find('MODEL        1') != -1:
                break
            else:
                current_header += line
        file_hand.close()
        return current_header
    except OSError as err:
        print("OS error: {0}".format(err))
    except ValueError:
        print("The file is damaged")

def merge_the_header(file_pdb, string_header, fileout):
    """
    This function merge the string header extracted from the original
    pdb file to the superimposed new pdb file. This way we restore the
    secondary structure visualitzation.
    """
    try:
        file_in = open(file_pdb, "r")
    except OSError as err:
        print("OS error: {0}".format(err))
    file_out = open(fileout, "w")
    file_out.write(string_header)
    for line in file_in:
        file_out.write(line)

def pdb_code_check(code):
    """Checks if the code provided have 4 letters"""
    try:
        return len(code) == 4
    except TypeError:
        return False
    except AttributeError:
        return False




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
    if pdb_code_check(code):
        pdbl = pdb.PDBList()
        parser = pdb.PDBParser(QUIET=True)
        try:
            structure = parser.get_structure(code, pdbl.retrieve_pdb_file(code, \
                pdir="pdbfiles/"))
        except urllib.error.URLError:
            sys.stderr.write("There is no a structure with the pdb code {} \
in the database\n".format(code))
        if nmr_check(structure):
            return structure
        else:
            return ValueError
    else:
        return ValueError
