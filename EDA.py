import os
import argparse as arg
import module_david as mdl
import Bio.PDB as pdb
import project_jf as calc
from scipy import linalg

parse_args = arg.ArgumentParser(description="This program performs the \
                                essential dynamics(ED) analysis of a protein, \
                                given a NMR solved structure or an MD \
                                trajectory. It provides two \
                                interfaces, a graphical and a command line \
                                the default is the graphical, and the command \
                                line can be used specifiyin the -ng option")
parse_args.add_argument('-ng', dest="graphical", action="store_true",
                        default=False, help="Use the command line interface, \
                        if it is not specified the graphical interface will \
                        be called ")
parse_args.add_argument('-i', '--input', dest="infile", action="store",
                        default=None, help="Input NMR pdb file or code or MD \
                        trajectories")
parse_args.add_argument('-m', '--mode', dest='mode', action="store",
                        choices=['NMR', 'MD'], help="Select wether you want to \
                        perform the ED analysis on NMR or MD data")
parse_args.add_argument('-a', '--atoms', dest='atom', action="store",
                        choices=['CA', 'Back', 'all'], help="Select wich atoms \
                        you want to perform the ED analysis, CA means only \
                        Carbon alpha atoms, Back means backbone atoms(CA, N, \
                        C, O) or all, this later option is discouraged since \
                        it may be computationally very expensive and take more \
                        time")
parse_args.add_argument('-v', '-verbose', dest='verb', action="store_true",
                        default=False, help="Verbose description of program \
                        actions")
options = parse_args.parse_args()

pathname = 'pdbfiles/'
pathplots = 'plots/'
if not os.path.exists(pathname):
    os.mkdir(pathname)
if not os.path.exists(pathplots):
    os.mkdir(pathplots)

if not options.graphical:
    import interface
    interface.quit_window()

if mdl.pdb_code_check(options.infile):
    pdb_id = options.infile
    pdbfile = pdb_id + '.ent'
else:
    pdbfile = options.infile
if not os.path.exists(pathname+pdbfile):
    pdbobj = pdb.PDBList()
    pdbfile = pdbobj.retrieve_pdb_file(pdb_id, pdir=pathname)
pdbalignedfile = pdb_id + 'align.pdb'
pdb_superimp = pathname + pdb_id + 'superimp.pdb'
parser = pdb.PDBParser(QUIET=options.verb)
structure = parser.get_structure(pdb_id, pdbfile)

atom_list = []
if options.atom == 'CA':
    atom_list = ['CA']
elif options.atom == 'Back':
    atom_list = ['N', 'CA', 'C', 'O']

N = 0
# Calculate the number of residues that are aa
for residue in structure[0].get_residues():
    if residue.has_id('CA'):
        N += 1
# Calculate the number of atoms to study, this will depend on which kind
# of atom the users wants to focus on(CA, backbone or all)
if atom_list != []:
    N *= len(atom_list)
# Calculate the number of configurations available, in NMR this is the
# number of models, in MD the number of time instants
n = len(structure)

head = mdl.store_header_text(pdbfile)
structure = calc.superimpose_models(structure, atom_list)
io = pdb.PDBIO()
io.set_structure(structure)
io.save(pdb_superimp)
mdl.merge_the_header(pdb_superimp, head, pathname+pdbalignedfile)
os.remove(pdb_superimp)

if options.verb:
    print("Calculating means and coordinates")
(array_stored, means) = calc.createcordsarray(structure, N, atom_list)
if options.verb:
    print("Calculating covariance matrix")
C = calc.cal_cov(array_stored, means)
if options.verb:
    print("Calculating eigenvalues and eigenvectors")
evl, evc = linalg.eigh(C)
if options.verb:
    print("Plotting eigenvalues")
plot = calc.plot_eig(evl, n, pathplots, pdb_id)
