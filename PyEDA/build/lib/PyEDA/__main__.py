#!/usr/bin/python3

import os
import sys
import argparse as arg
import Bio.PDB as pdb
from . import helper_module as mdl
from . import edanalysis as eda


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
parse_args.add_argument('-c', '--code', dest="code", action="store",
                        default=None, help="Input pdb code for the protein you \
                        want to study")
parse_args.add_argument('-i', '--input', dest="infile", action="store",
                        default=None,
                        help="Input pdb file with either NMR or MD data")
parse_args.add_argument('-m', '--mode', dest='mode', action="store",
                        choices=['NMR', 'MD'], help="Select wether you want to \
                        perform the ED analysis on NMR or MD data",
                        default=None)
parse_args.add_argument('-a', '--atoms', dest='atom', default='CA',
                        action="store", choices=['CA', 'Back', 'all'],
                        help="Select wich atoms you want to perform the \
                        ED analysis, CA means only Carbon alpha atoms,  \
                        Back means backbone atoms(CA, N, C, O) or all,  \
                        this later option is discouraged since it may be \
                        computationally very expensive and take more \
                        time")
parse_args.add_argument('-d', '--directory', dest='path', action='store',
                        help="Path to the directory where you want to work, if \
                        not specified, current working directory will be used")
parse_args.add_argument('-v', '-verbose', dest='verb', action="store_true",
                        default=False, help="Verbose description of program \
                        actions")
parse_args.add_argument('--eigvl-num', dest='eigvl', action='store', default=20,
                        help="Number of eigenvalues to plot (default is 20)")
parse_args.add_argument('--time-max', dest='time', action='store', default=1,
                        help="Maximum time during which the trajectories of a \
                        certain eigenvector")
parse_args.add_argument('--time-step', dest='step', action='store',
                        default=0.1, help="Time step for the generated \
                        trajectories")
parse_args.add_argument('--eigvc-num', dest='eigvc', action='store', default=1,
                        help="Number of eigenvector for which new trajectories \
                        will be generated, default is 1 since it is a \
                        computationally expensive step")
options = parse_args.parse_args()

path = ''
if options.path:
    path = options.path
pathname = path + 'pdbfiles/'
pathplots = path + 'plots/'

if not os.path.exists(pathname):
    os.mkdir(pathname)
if not os.path.exists(pathplots):
    os.mkdir(pathplots)

if not options.graphical:
    from . import interface as inter
    interface = inter.EDA_app()
    interface.mainloop()
    sys.exit("The application has been closed.\n")

if options.mode is None:
    raise ValueError('Please specify a mode with the -m option')
if options.code is None:
    raise ValueError('Please specify a PDB code for your protein with the '
                     "-c option")
if options.mode == 'MD' and options.infile is None:
    raise ValueError(('To use the MD mode you need to input a trajectory in '
                      'a PDB file'))
if mdl.pdb_code_check(options.code):
    pdb_id = options.code
else:
    raise ValueError('Input code is not a PDB code')

if options.infile:
    pdbfile = options.infile
else:
    pdbfile = 'pdb'+pdb_id + '.ent'

pdbalignedfile = pdb_id + 'align.pdb'
pdb_superimp = pathname + pdb_id + 'superimp.pdb'

if not os.path.exists(pathname+pdbfile):
    pdbobj = pdb.PDBList()
    pdbobj.retrieve_pdb_file(pdb_id, pdir=pathname)

if not (pdbfile.endswith('pdb') or pdbfile.endswith('ent')):
    raise ValueError(('Your input file is not a valid PDB file, please use a '
                      'pdb or ent file'))

atom_list = []
if options.atom == 'CA':
    atom_list = ['CA']
elif options.atom == 'Back':
    atom_list = ['N', 'CA', 'C', 'O']

if options.verb:
    print("Initializing analysis information")

if options.mode == 'MD':
    pdbref = pdb.PDBList()
    ref_file = pdbref.retrieve_pdb_file(pdb_id, pdir=pathname)
    parser = pdb.PDBParser(QUIET=True)
    reference = parser.get_structure(pdb_id+'ref', ref_file)
    head = mdl.store_header_text(pathname+pdbfile)
    ED = eda.EDAnalysis(pdb_id, options.mode, atom_list, pathname+pdbfile,
                        reference=reference)
else:
    ED = eda.EDAnalysis(pdb_id, options.mode, atom_list, pathname+pdbfile)

if options.eigvl < 1 or options.eigvl > ED.N:
    raise ValueError('Number of eigenvalues to plot must be between 1 and N')
if options.eigvc < 1 or options.eigvc > ED.N:
    raise ValueError('Number of eigenvectors must be between 1 and N')
if options.time < 0:
    raise ValueError('Final time for the generated trajectories must be > 0')

if options.verb:
    print("Superimposing structures to a reference")
ED.superimpose_models()

if options.mode == 'NMR':
    head = mdl.store_header_text(pathname+pdbfile)

io = pdb.PDBIO()
io.set_structure(ED.structure)
io.save(pdb_superimp)
mdl.merge_the_header(pdb_superimp, head, pathname+pdbalignedfile)
os.remove(pdb_superimp)

if options.verb:
    print("Calculating means and coordinates")
ED.createcordsarray()

if options.verb:
    print("Calculating covariance matrix")
    print("Calculating eigenvalues and eigenvectors")
ED.cal_cov()

if options.verb:
    print("Plotting eigenvalues")
n_plot = options.eigvl
if ED.n < n_plot:
    n_plot = ED.n
plot = ED.plot_eig(n_plot, pathplots)

if options.verb:
    print("Generating eigenvector trajectories")
moved = ED.move_structure(options.time, options.eigvc, pathname,
                          step=options.step)
new_moved = moved[:-5]+'.pdb'
mdl.merge_the_header(moved, head, new_moved)
os.remove(moved)

if options.atom != 'all':
    if options.verb:
        print("Generating RMSD plot for eigenvectors")
    plot = ED.RMSD_res_plot(options.eigvc, pathplots)
