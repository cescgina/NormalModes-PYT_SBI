
import Bio.PDB as pdb
import os
import numpy as np

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

structure['header']
io = pdb.PDBIO()
io.set_structure(structure)
io.save(pathname+pdbalignedfile)
