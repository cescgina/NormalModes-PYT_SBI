import copy
import Bio.PDB as pdb
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy import linalg
matplotlib.use("TkAgg")


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

    def plot_eig_wosv(self, n_plot, fig=None):
        """ """
        if fig is None:
            fig = plt.figure()
        valid_evl = self.eigvl[-1:-n_plot-1:-1] / 100
        plt.plot(range(1, n_plot+1), valid_evl)
        plt.xlabel('Eigenvector index')
        plt.ylabel('Eigenvalue ($nm^2$)')
        # plt.axis([0, n, 0, 3])
        return fig

    def is_NMR_struct(self):
        """ Function to ensure the loaded structure is a NMR"""
        return 'nmr' in self.structure.header['structure_method'].lower()

    def check_mode(self):
        """ """
        if self.__mode is None:
            raise WrongModeException('empty')
        if self.__mode == 'NMR' and not self.is_NMR_struct():
            raise WrongModeException(self.__mode)
        elif self.__mode == 'MD' and self.is_NMR_struct():
            raise WrongModeException(self.__mode)

    def move_structure(self, t_max, evc, pathname, t_min=0, step=0.1):
        """ Returns a new structure with the coordinates calculated from a
        certain eigenvector"""
        structure_moved = copy.deepcopy(self.structure)
        filename = ''.join([pathname, self.__PDBid, '_evc', str(evc), '_0.pdb'])
        io = pdb.PDBIO()
        io.set_structure(structure_moved)
        io.save(filename, OneChainSelect())
        pcord = np.dot(self.eigvc[:, -evc], (self.coords_array[0] -
                                             self.means[0, :]))
        image_list = [filename]
        nsteps = int((t_max - t_min) / step)
        for t in np.linspace(t_min, t_max, num=nsteps):
            if int(t/step) == 0:
                # This block maybe could be avoided and use the means as
                # reference
                continue
            eig_move = t * pcord * self.eigvc[:, -evc] + self.means
            j = 0
            for residue in structure_moved[0].get_residues():
                if residue.has_id('CA'):
                    for atoms in self.__atom:
                        residue[atoms].set_coord(eig_move[0][j:j+3])
                        j += 3
            filename = ''.join([pathname, self.__PDBid, '_evc', str(evc), '_',
                                str(int(t/step)), '.pdb'])
            io = pdb.PDBIO()
            io.set_structure(structure_moved)
            io.save(filename, OneChainSelect())
            image_list.append(filename)
        return image_list

    def RMSD_res_plot(self, evc, pathplots, fig=None):
        """ """
        if fig is None:
            fig = plt.figure()
        for evcn in range(1, evc+1):
            pcord = 0
            pmin = 10000000
            pmax = -1000000
            for ind in range(self.n):
                # Look for the maximum and minimum translation
                pcord = np.dot(self.eigvc[:, -evcn], (self.coords_array[ind] -
                                                      self.means[0, :]))
                if pcord > pmax:
                    pmax = pcord
                elif pcord < pmin:
                    pmin = pcord
            eig_move_max = pmax * self.eigvc[:, -evcn] + self.means
            eig_move_min = pmin * self.eigvc[:, -evcn] + self.means
            pcord = np.dot(self.eigvc[:, -evcn], (self.coords_array[0] -
                                                  self.means[0, :]))
            # eig_move = pcord * self.eigvc[:, -evcn] + self.means
            step = len(self.__atom)*3
            nres = int(self.N / step) - 15 + 1 + 1
            RMSD_list = np.zeros(nres)
            j = 0
            i = 0

            for residue in self.structure[0].get_residues():
                if residue.has_id('CA'):
                    j_final = j+7*step
                    j_init = j-7*step
                    if j_final > self.N:
                        break
                    elif j_init < 0:
                        j += step
                        continue
                    else:
                        # RMSDvl = eig_move[0][j_init:j_final] - \
                        #      self.coords_array[0][j_init:j_final]
                        RMSDvl = eig_move_max[0][j_init:j_final] - \
                           eig_move_min[0][j_init:j_final]
                        j += step
                        RMSD = np.sqrt(np.sum(RMSDvl**2)/len(RMSDvl))
                        RMSD_list[i] = RMSD
                        i += 1
            plt.plot(range(7, nres+7), RMSD_list,
                     label="EV {:d}".format(evcn))
        plt.ylabel('RMSD ($\AA$)')
        plt.xlabel('Residue number')
        # plt.axis([0, n, 0, 3])
        plt.legend(loc='best', frameon=False)
        filename = ''.join([pathplots, 'eig_', str(evc), '_', self.__PDBid,
                            '_resplot.png'])
        fig.savefig(filename, bbox_inches='tight', dpi=300)
        return fig


class OneChainSelect(pdb.Select):
    """
    Custom class derived from Bio.PDB.Select to write only one model to a PDB
    file, used only to visualize the trajectories of the eigenvectors
    """
    def accept_model(self, model):
        return model.get_id() == 0


class WrongModeException(Exception):
    def __init__(self, mode):
        self.mode = mode

    def __str__(self):
        """ 2"""
        return ("Your mode selection is {:s}, please input an appropiate "
                "structure").format(self.mode)
