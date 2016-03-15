#!/usr/bin/python3

import os
import sys
import argparse as arg
import module_david as mdl
import Bio.PDB as pdb
import edanalysis as eda
# import __main__
# __main__.pymol_argv = ['pymol', '-qc']
# import pymol
# pymol.finish_launching()


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
                        perform the ED analysis on NMR or MD data")
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
    import interface as inter
    interface = inter.EDA_app()
    interface.mainloop()
    sys.exit("The application has been closed.\n")

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

ED = eda.EDAnalysis(pdb_id, options.mode, atom_list, pathname+pdbfile)

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
n_plot = 30
if ED.n < n_plot:
    n_plot = ED.n
plot = ED.plot_eig(n_plot, pathplots)

# if options.verb:
#     print("Generating eigenvector trajectories")
# image_list = ED.move_structure(10, 1, pathname)

if options.verb:
    print("Generating RMSD plot for eigenvectors")
plot = ED.RMSD_res_plot(4, pathplots)
# for file_img in image_list:
#     sname = file_img.rstrip(".pdb")
#     pymol.cmd.load(file_img, sname)
#     pymol.cmd.disable("all")
#     pymol.cmd.enable(sname)
#     pymol.cmd.show("cartoon")
#     pymol.cmd.png(sname+".png")
# # Get out!
# pymol.cmd.quit()
