"""This module is part of the EDA(essential dynamics analysis) program, it
provides several classes to perform the EDA.

Classes:
    EDAnalysis -> It is the main class, provides an easy interface to
    carry out ED
    WrongModeException -> Specific Exception class that indicates a misuse
    of the EDA program
"""
import copy
import Bio.PDB as pdb
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from scipy import linalg


class EDAnalysis:
    """Implements an essential dynamics analysis(EDA).

    Methods:

        __init__ -> Set the basic attributes of the EDAnalysis
        get_id -> Return the id of the object (as PDB code)
        get_mode -> Return the mode of the object
        get_atom -> Return the atom list chosen for the analysis
        superimpose_models- > Superimpose two or more structures
        createcordsarray -> Calculate the coordinates for the diferent models
        cal_cov -> Calculate the covariance matrix
        plot_eig -> Plot the values of the principal eigenvectors of the
             covariance matrix
        plot_eig_wosv -> Same as plot_eig but adapted to some requirements of
            the Tkinter GUI
        is_NMR_struct -> Check if the structure loaded is was obtained by NMR
        check_mode -> Check if the mode selected is consisted with the data
        move_structure -> Move the structure along a certain eigenvector
        RMSD_res_plot -> Create a plot with the distance of the residues along
            some eigenvectors

    Attributes:

        PDBid(private) -> String with the PDB code of the structure anlyzed
        mode(private) -> String with the mode of the anlysis (for now NMR or MD)
        atom(private) -> List with the atoms to study
        structure -> Bio.PDB.Structure object on which the analysis is performed
        reference -> Bio.PDB.Structure object with a X-Ray structure
                that will serve as a reference model, this should only be used
                with the MD mode
        n -> Number of models or trajectories (in NMR or MD modes,respectively)
        N -> Number of coordinates of the analysis
        coords_array -> nxN Numpy array with the N coordinates for the n models
        means -> 1xN Numpy array with the means of the N coordinates over the n
            models
        C -> NxN Numoy array that contains the covariance matrix
        eigvc -> Numpy array with the eigenvectors of the covariance matrix
        eigvl -> Numpy array with the eigenvalues of the covariance matrix
    """

    def __init__(self, PDBid, mode, atom, pdbfile, reference=None):
        """Set the basic attributes of the EDAnalysis class.

        Set the basic attributes of the EDAnalysis class, calculate the number
        of coordinates of interest, the number of models, load the PDB file
        with the structure to analyze and check the consistency of the data
        provided by the user

        Args:
            PDBid is a string with a PDB code
            mode is a string with value NMR or MD
            atom is a list with the particular atoms to focus on the analysis
            pdbfile is the path to pdbfile to load
        """

        self.__PDBid = PDBid
        self.__mode = mode
        self.__atom = atom
        if reference is not None:
            if self.__mode != 'MD':
                raise WrongModeException(self.__mode)
        self.reference = reference

        parser = pdb.PDBParser(QUIET=True)
        self.structure = parser.get_structure(self.__PDBid, pdbfile)
        N = 0
        # Calculate the number of residues that are aa
        for residue in self.structure[0].get_residues():
            if residue.has_id('CA'):
                if self.__atom == []:
                    N += len(list(residue.get_atom()))
                else:
                    N += 1
        # Calculate the number of atoms to study, this will depend on which
        # atoms the users wants to focus(CA, backbone or all)
        if self.__atom != []:
            N *= len(self.__atom)*3
        else:
            N *= 3
        # Calculate the number of configurations available, in NMR this is the
        # number of models, in MD the number of time instants
        self.n = len(self.structure)
        self.N = N
        self.check_mode()

    def get_id(self):
        """Return the id of the EDAnalysis as a string."""
        return self.__PDBid

    def get_mode(self):
        """Return the selected mode of the EDAnalysis as a string."""
        return self.__mode

    def get_atom(self):
        """Return the list of atoms of interest to the EDAnalysis."""
        return self.__atom

    def superimpose_models(self, reference_model=0):
        """Superimpose two or more structures.

        Superimpose two or more structures by using the Bio.PDB.Superimposer
        class.

        Args:

            reference_model is an int with the reference model (default=0)
        """

        ref_model = self.structure[reference_model]
        if self.reference is not None:
            ref_model = self.reference[0]

        for alt_model in self.structure:
            ref_atoms = []
            alt_atoms = []
            # Iterate over the structure method to obtain all atoms of interest
            # for the analysis and superimposes the structures using them
            for (ref_chain, alt_chain) in zip(ref_model, alt_model):
                for ref_res, alt_res in \
                        zip(ref_chain, alt_chain):
                    # assert ref_res.resname == alt_res.resname, \
                        # "{:s} is not equal to {:s}".format(ref_res.resname,
                                                            # alt_res.resname)
                    # assert ref_res.id == alt_res.id, \
                        # "{:s} is not equal to {:s}".format(ref_res.id,
                                                            # alt_res.id)
                    # CA = alpha carbon
                    if ref_res.has_id('CA'):
                        if self.__atom == []:
                            ref_atoms.extend(list(ref_res.get_atom()))
                            alt_atoms.extend(list(alt_res.get_atom()))
                        else:
                            for atoms in self.__atom:
                                try:
                                    ref_atoms.append(ref_res[atoms])
                                    alt_atoms.append(alt_res[atoms])
                                except KeyError:
                                    raise KeyError(('Your input data is '
                                                    'misssing information for '
                                                    '{:s} atoms. Input more '
                                                    'complete data or select a'
                                                    ' smaller set of '
                                                    'atoms').format(atoms))

            # Align these paired atom lists:
            super_imposer = pdb.Superimposer()
            super_imposer.set_atoms(ref_atoms, alt_atoms)

            if ref_model.get_full_id() == alt_model.get_full_id():
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
        """Calculate the coordinates for the diferent models."""
        array_stored = np.zeros((self.n, self.N))
        means = np.zeros((1, self.N))
        i = 0
        for model in self.structure:
            j = 0
            for residue in model.get_residues():
                if residue.has_id('CA'):
                    if self.__atom != []:
                        for atoms in self.__atom:
                            try:
                                array_stored[i][j:j+3] = residue[atoms].get_coord()
                                means[0][j:j+3] += residue[atoms].get_coord()
                                j += 3
                            except KeyError:
                                raise KeyError(('Your input data is misssing '
                                                'information for {:s} atoms'
                                                '. Input more complete data '
                                                'or select a smaller set of '
                                                'atoms').format(atoms))
                    else:
                        for atoms in residue.get_atom():
                                array_stored[i][j:j+3] = residue[atoms.id].get_coord()
                                means[0][j:j+3] += residue[atoms.id].get_coord()
                                j += 3
            i += 1
        means *= (1/self.n)
        self.coords_array = array_stored
        self.means = means

    def cal_cov(self):
        """Calculate the covariance matrix.

        Calculate the covariance matrix from the coordinates caculated with the
        createcordsarray method, the diagonalize the covariance matrix using
        scipy.linalg.eigh
        """

        C = np.zeros((self.N, self.N))
        for x in self.coords_array:
            x_prod = x - self.means
            C += np.outer(x_prod, x_prod)
        C *= (1/self.n)
        self.C = C
        self.eigvl, self.eigvc = linalg.eigh(C)

    def plot_eig(self, n_plot, pathplots, fig=None):
        """Plot the values of the principal eigenvectors of the cov. matrix.

        Args:
            n_plot -> Int with the number of eigvalues to plot
            pathplots -> Directory where the plot will be stored
            fig -> Matplotlib figure handle, the default in None, so a new
                figure handle will be created

        Returns the figure handle that contains the plot
        """

        if fig is None:
            fig = plt.figure()
        valid_evl = self.eigvl[-1:-n_plot-1:-1] / 100
        plt.plot(range(1, n_plot+1), valid_evl)
        plt.xlabel('Eigenvector index')
        plt.ylabel('Eigenvalue ($nm^2$)')
        fig.savefig(pathplots+'eig_'+self.__PDBid+'_plot.png',
                    bbox_inches='tight', dpi=300)
        return fig

    def plot_eig_wosv(self, n_plot, fig=None):
        """Essentially the same as plot_eig, just adapted to the Tkinter GUI."""
        if fig is None:
            fig = plt.figure()
        valid_evl = self.eigvl[-1:-n_plot-1:-1] / 100
        plt.plot(range(1, n_plot+1), valid_evl)
        plt.xlabel('Eigenvector index')
        plt.ylabel('Eigenvalue ($nm^2$)')
        return fig

    def is_NMR_struct(self):
        """Check wether the loaded structure is a NMR"""
        return 'nmr' in self.structure.header['structure_method'].lower()

    def check_mode(self):
        """Check the consistency betweent the mode chosen by the user and the
        data provided.
        """

        if self.__mode is None:
            raise WrongModeException('empty')
        if self.__mode == 'NMR' and not self.is_NMR_struct():
            raise WrongModeException(self.__mode)
        elif self.__mode == 'MD' and self.is_NMR_struct():
            raise WrongModeException(self.__mode)

    def move_structure(self, t_max, evc, pathname, t_min=0, step=0.1):
        """Move the structure along a certain eigenvector.

        Project the coordinates of the stucture onto an eigenvector and then
        uses this projection to move the structure along this eigenvector,
        creating a PDB file with the now moved coordinates

        Args:
            t_max -> Int with the final time of the trajectory
            evc -> Int with the eigenvector to use, raise a ValueError if it is
                not between 1 and N
            pathname -> Directory where to store the pdb files generated
            t_min -> Int with initial time of the trajectory (default=0)
            step -> Float with the time-step of the trajectory (default=0.1)

        Return the name of the file containing the new coordinates generated
        as diferent models
        """

        if evc < 1 or evc > self.N:
            raise ValueError('Eigenvector index has to be between 1 and N')
        structure_moved = copy.deepcopy(self.structure)
        pcord = np.dot(self.eigvc[:, -evc], (self.coords_array[0] -
                                             self.means[0, :]))
        nsteps = int((t_max - t_min) / step)
        for i in range(self.n):
            structure_moved.detach_child(i)
        idnum = 0
        for t in np.linspace(t_min, t_max, num=nsteps):
            new_model = copy.deepcopy(self.structure[0])
            new_model.full_id = (self.structure.id, idnum)
            new_model.id = idnum
            new_model.serial_num = idnum+1
            idnum += 1
            eig_move = t * pcord * self.eigvc[:, -evc] + self.means
            j = 0
            for residue in new_model.get_residues():
                if residue.has_id('CA'):
                    if self.__atom != []:
                        for atoms in self.__atom:
                            try:
                                residue[atoms].set_coord(eig_move[0][j:j+3])
                                j += 3
                            except KeyError:
                                raise KeyError(('Your input data is misssing '
                                                'information for {:s} atoms'
                                                '. Input more complete data '
                                                'or select a smaller set of '
                                                'atoms').format(atoms))
                    else:
                        for atoms in residue.get_atom():
                            residue[atoms.id].set_coord(eig_move[0][j:j+3])
                            j += 3
            structure_moved.add(new_model)
        filename = ''.join([pathname, self.__PDBid, '_traj_evc', str(evc),
                            'i.pdb'])
        self.structure_moved = structure_moved
        io = pdb.PDBIO()
        io.set_structure(structure_moved)
        # io.save(filename, ChainSelect(nsteps))
        io.save(filename)
        return filename

    def RMSD_res_plot(self, evc, pathplots, fig=None, origin=None):
        """Create a plot with the distance of the residues along
        some eigenvectors

        Project the structure over some eigenvectors and calculte new
        coordinates, calculate the RMSD between the new and old coordinates by
        using a shifting-window RMS method (with a window of 15 residues) and
        plot the RMSD against the residue number

        Args:
            evc -> Int with the maximum eigenvector to include in the plot
            pathplots -> Directory to store the plot
            fig -> Matplotlib figure handle, the default in None, so a new
                figure handle will be created

        Returns the figure handle that contains the plot
        """
        if fig is None:
            fig = plt.figure()

        if evc < 1 or evc > self.N:
            raise ValueError('Eigenvector index has to be between 1 and N')

        step = 1
        if self.__atom != []:
            step = len(self.__atom)*3
        nres = int(self.N / step) - 15 + 1 + 1

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
            RMSD_list = np.zeros(nres)
            j = 0
            i = 0

            for residue in self.structure[0].get_residues():
                # Calculate the RMSD for each residue using the shifting-window
                # RMS method with a window size of 15
                if residue.has_id('CA'):
                    j_final = j+7*step
                    j_init = j-7*step
                    if j_final > self.N:
                        break
                    elif j_init < 0:
                        j += step
                        continue
                    else:
                        RMSDvl = eig_move_max[0][j_init:j_final] - \
                           eig_move_min[0][j_init:j_final]
                        j += step
                        RMSD = np.sqrt(np.sum(RMSDvl**2)/int(len(RMSDvl)/step))
                        RMSD_list[i] = RMSD
                        i += 1
            plt.plot(range(7, nres+7), RMSD_list,
                     label="EV {:d}".format(evcn))
        plt.ylabel('RMSD ($\AA$)')
        plt.xlabel('Residue number')
        plt.legend(loc='best', frameon=False)
        filename = ''.join([pathplots, 'eig_', str(evc), '_', self.__PDBid,
                            '_resplot.png'])
        if origin != 'interface':
            fig.savefig(filename, bbox_inches='tight', dpi=300)
        return fig


class WrongModeException(Exception):
    """Specific Exception class that indicates a misuse of the EDA program."""

    def __init__(self, mode):
        self.mode = mode

    def __str__(self):
        return ("Your mode selection is {:s}, please input an appropiate "
                "structure").format(self.mode)
