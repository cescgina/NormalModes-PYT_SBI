
import Bio.PDB as pdb
import os
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import module_david as mdl


class WrongModeException(Exception):
    def __init__(self, input_class, mode):
        self.input_class = input_class
        self.mode = mode

    def __str__(self):
        """ 2"""
        return "You have selected the wrong mode, your input is %s \
                and you have selected mode %s" % (self.input_class, self.mode)


def is_NMR_struct(structure):
    """ Function to ensure the loaded structure is a NMR"""
    return 'nmr' in structure.header['structure_method'].lower()


def superimpose_models(structure, atom_list):
    """ Function to sumperimpose all the models of the NMR pdb file into
    a reference model, which is assumed to be the first model of the ensemble
    """
    print("Superimposing different models onto the first one")
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
                    if atom_list == []:
                        ref_atoms.extend(list(ref_res.get_atoms()))
                        alt_atoms.extend(list(alt_res.get_atoms()))
                    else:
                        for atom in atom_list:
                            ref_atoms.append(ref_res[atom])
                            alt_atoms.append(alt_res[atom])

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
        print("RMS(Refernce model, model %i) = %0.2f" %
              (alt_model.id, super_imposer.rms))
    return structure


def createcordsarray(structure, N, atom_list):
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
                for atom in atom_list:
                    array_stored[i][j:j+3] = residue[atom].get_coord()
                    means[0][j:j+3] += residue[atom].get_coord()
                    j += 3
    means *= (1/n)
    return (array_stored, means)


def cal_cov(array_stor, means):
    """ Helper function to calculate the covariance matrix"""
    (n, N) = array_stor.shape
    C = np.zeros((N, N))
    for k in range(n):
        for ii in range(N):
            for ij in range(ii, N):
                C[ii][ij] += (1 / n) * (array_stor[k][ii] - means[0][ii]) * \
                                (array_stor[k][ij] - means[0][ij])
                C[ij][ii] = C[ii][ij]
    return C

# pdb_id = '1jm7'
pdb_id = '1msf'
pdbfile = pdb_id + '.ent'
pdbalignedfile = pdb_id + 'align.pdb'
pathname = 'pdbfiles/'
pdb_superimp = pathname + pdb_id + 'superimp.pdb'
if not os.path.exists('pdbfiles'):
    os.mkdir('pdbfiles')
if not os.path.exists(pathname+pdbfile):
    pdbobj = pdb.PDBList()
    pdbfile = pdbobj.retrieve_pdb_file(pdb_id, pdir=pathname)

# atom_list = ['CA']
atom_list = ['N', 'CA', 'C', 'O']
parser = pdb.PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdbfile)
header = parser.get_header()
if not is_NMR_struct(structure):
    raise WrongModeException(structure.header['structure_method'], 'NMR')
N = 0
for residue in structure[0].get_residues():
    if residue.has_id('CA'):
        N += 1

if atom_list != []:
    N *= len(atom_list)
n = len(structure)

head = mdl.store_header_text(pdbfile)
structure = superimpose_models(structure, atom_list)
io = pdb.PDBIO()
io.set_structure(structure)
io.save(pdb_superimp)
mdl.merge_the_header(pdb_superimp, head, pathname+pdbalignedfile)

print("Calculating means and coordinates")
(array_stored, means) = createcordsarray(structure, N, atom_list)
print("Calculating covariance matrix")
C = cal_cov(array_stored, means)
print("Calculating eigenvalues and eigenvectors")
evl, evc = linalg.eigh(C)

valid_evl = evl[-1:-n-1:-1] / 100
plt.plot(range(1, n+1), valid_evl)
plt.xlabel('Eigenvector index')
plt.ylabel('Eigenvalue ($nm^2$)')
plt.axis([0, n, 0, 3])
plt.savefig('eig_plot.png', bbox_inches='tight')
plt.show()
