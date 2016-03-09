
import Bio.PDB as pdb
import os
import numpy as np
from scipy import linalg

# pdb_id = '1jm7'
pdb_id = '1joy'
pdbfile = pdb_id + '.ent'
pdbalignedfile = pdb_id + 'align.pdb'
pathname = 'pdbfiles/'
if not os.path.exists('pdbfiles'):
    os.mkdir('pdbfiles')
if not os.path.exists(pathname+pdbfile):
    pdbobj = pdb.PDBList()
    pdbfile = pdbobj.retrieve_pdb_file(pdb_id, pdir=pathname)

parser = pdb.PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdbfile)
header = parser.get_header()

N = 0
for residue in structure[0].get_residues():
    if residue.has_id('CA'):
        N += 1


def Superimpose_models(structure):
    """ Function to sumperimpose all the models of the NMR pdb file into
    a reference model, which is assumed to be the first model of the ensemble
        Writes a new pdb with the superimposed coordinates
    """
    ref_model = structure[0]
    for alt_model in structure:
        # Code from
        # http://www.warwick.ac.uk/go/peter_cock/python/protein_superposition/
        # Build paired lists of c-alpha atoms, ref_atoms and alt_atoms
        ref_atoms = []
        alt_atoms = []
        for (ref_chain, alt_chain) in zip(ref_model, alt_model):
            for ref_res, alt_res in \
                    zip(ref_chain, alt_chain):
                assert ref_res.resname == alt_res.resname
                assert ref_res.id == alt_res.id
                # CA = alpha carbon
                if ref_res.has_id('CA'):
                    ref_atoms.append(ref_res['CA'])
                    alt_atoms.append(alt_res['CA'])

        # Align these paired atom lists:
        super_imposer = pdb.Superimposer()
        super_imposer.set_atoms(ref_atoms, alt_atoms)

        if ref_model.id == alt_model.id:
            # Check for self/self get zero RMS, zero translation
            # and identity matrix for the rotation.
            assert np.abs(super_imposer.rms) < 0.0000001
            assert np.max(np.abs(super_imposer.rotran[1])) < 0.000001
            assert np.max(np.abs(super_imposer.rotran[0]) - np.identity(3)) \
                < 0.000001
        else:
            # Update the structure by moving all the atoms in
            # this model (not just the ones used for the alignment)
            super_imposer.apply(alt_model.get_atoms())
        print("RMS(first model, model %i) = %0.2f" %
              (alt_model.id, super_imposer.rms))
    # Due to the way Bio.PDB.PDBIO is written, the pdbalignedfile
    # does not contain any header or trailer information, maybet it is not
    # necessary to write the modified structure at all, since the
    # superimposition changes are applied in the Structure object
    io = pdb.PDBIO()
    io.set_structure(structure)
    io.save(pathname+pdbalignedfile)


def createcordsarray(structure, N):
    """ Helper function to extract the coordinates of the atoms to calculate
    the covariance matrix afterwards
    """
    n = len(structure)
    array_stored = np.zeros((n, 3*N))
    means = np.zeros((1, 3*N))
    for i in range(n):
        j = 0
        for residue in structure[i].get_residues():
            if residue.has_id('CA'):
                array_stored[i][j:j+3] = residue['CA'].get_coord()
                means[0][j:j+3] += residue['CA'].get_coord()
                j += 3
    means *= (1/n)
    return (array_stored, means)


def cal_cov(array_stor, means, structure):
    """ Helper function to calculate the covariance matrix"""
    N = means.size
    C = np.zeros((N, N))
    n = len(structure[0])
    for k in range(n):
        for ii in range(N):
            for ij in range(ii, N):
                C[ii][ij] += array_stored[k][ii]*array_stor[k][ij] - \
                    array_stor[k][ii]*means[0][ij]-array_stor[k][ij] * \
                    means[0][ii] + means[0][ij]*means[0][ii]
                C[ii][ij] *= (1/n)
                C[ij][ii] = C[ii][ij]
    return C

(array_stored, means) = createcordsarray(structure, N)
C = cal_cov(array_stored, means, structure)
evl, evc = linalg.eigh(C)
