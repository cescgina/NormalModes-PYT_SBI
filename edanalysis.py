import Bio.PDB as pdb
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


class EDAnalysis:
    """Implements an essential dynamics analysis..."""
    def __init__(self, PDBid, mode, atom, pdbfile):
        self.__PDBid = PDBid
        self.__mode = mode
        self.__atom = atom
        parser = pdb.PDBParser(QUIET=True)
        self.structure = parser.get_structure(self.__PDBid, pdbfile)
        N = 0
        # Calculate the number of residues that are aa
        for residue in self.structure[0].get_residues():
            if residue.has_id('CA'):
                N += 1
        # Calculate the number of atoms to study, this will depend on which
        # atoms the users wants to focus(CA, backbone or all)
        if self.__atom != []:
            N *= len(self.__atom)*3
        # Calculate the number of configurations available, in NMR this is the
        # number of models, in MD the number of time instants
        self.n = len(self.structure)
        self.N = N
        self.check_mode()

    def get_id(self):
        """ """
        return self.__PDBid

    def get_mode(self):
        """ """
        return self.__mode

    def get_atom(self):
        """ """
        return self.__atom

    def superimpose_models(self, reference=0):
        """ reference is an integer that indicates the reference model"""
        ref_model = self.structure[reference]
        for alt_model in self.structure:
            ref_atoms = []
            alt_atoms = []
            for (ref_chain, alt_chain) in zip(ref_model, alt_model):
                for ref_res, alt_res in \
                        zip(ref_chain, alt_chain):
                    assert ref_res.resname == alt_res.resname
                    assert ref_res.id == alt_res.id
                    # CA = alpha carbon
                    if ref_res.has_id('CA'):
                        if self.__atom == []:
                            ref_atoms.extend(list(ref_res.get_atoms()))
                            alt_atoms.extend(list(alt_res.get_atoms()))
                        else:
                            for atoms in self.__atom:
                                ref_atoms.append(ref_res[atoms])
                                alt_atoms.append(alt_res[atoms])

            # Align these paired atom lists:
            super_imposer = pdb.Superimposer()
            super_imposer.set_atoms(ref_atoms, alt_atoms)

            if ref_model.id == alt_model.id:
                # Check for self/self get zero RMS, zero translation
                # and identity matrix for the rotation.
                assert np.abs(super_imposer.rms) < 0.0000001
                assert np.max(np.abs(super_imposer.rotran[1])) < 0.000001
                assert np.max(np.abs(super_imposer.rotran[0]) -
                              np.identity(3)) < 0.000001
            else:
                # Update the structure by moving all the atoms in
                # this model (not just the ones used for the alignment)
                super_imposer.apply(alt_model.get_atoms())

    def createcordsarray(self):
        """ """
        array_stored = np.zeros((self.n, self.N))
        means = np.zeros((1, self.N))
        i = 0
        for model in self.structure:
            j = 0
            for residue in model.get_residues():
                if residue.has_id('CA'):
                    for atoms in self.__atom:
                        array_stored[i][j:j+3] = residue[atoms].get_coord()
                        means[0][j:j+3] += residue[atoms].get_coord()
                        j += 3
            i += 1
        means *= (1/self.n)
        self.coords_array = array_stored
        self.means = means

    def cal_cov(self):
        """ """
        C = np.zeros((self.N, self.N))
        for x in self.coords_array:
            x_prod = x - self.means
            C += np.outer(x_prod, x_prod)
        C *= (1/self.n)
        self.C = C
        self.eigvl, self.eigvc = linalg.eigh(C)

    def plot_eig(self, n_plot, pathplots, fig=None):
        """ """
        if fig is None:
            fig = plt.figure()
        valid_evl = self.eigvl[-1:-n_plot-1:-1] / 100
        plt.plot(range(1, n_plot+1), valid_evl)
        plt.xlabel('Eigenvector index')
        plt.ylabel('Eigenvalue ($nm^2$)')
        # plt.axis([0, n, 0, 3])
        fig.savefig(pathplots+'eig_'+self.__PDBid+'_plot.png',
                    bbox_inches='tight', dpi=300)
        return fig

    def is_NMR_struct(self):
        """ Function to ensure the loaded structure is a NMR"""
        return 'nmr' in self.structure.header['structure_method'].lower()

    def check_mode(self):
        """ """
        if self.__mode == 'NMR' and not self.is_NMR_struct():
            raise WrongModeException(self.__mode)
        elif self.__mode == 'MD' and self.is_NMR_struct():
            raise WrongModeException(self.__mode)


class WrongModeException(Exception):
    def __init__(self, input_class, mode):
        self.mode = mode

    def __str__(self):
        """ 2"""
        return "You have selected mode %s, please input an appropiate \
            structure " % (self.mode)
