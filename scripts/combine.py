from psfgen import PsfGen
from HyresBuilder import utils

RNA_topology, _ = utils.load_ff('RNA_m')
protein_topology, _ = utils.load_ff('protein_m')

gen = PsfGen()

gen.read_topology(RNA_topology)
gen.read_topology(protein_topology)

gen.add_segment(segid='PROA', pdbfile='chainA_cg.pdb')
gen.read_coords(segid='PROA', filename='chainA_cg.pdb')
gen.add_segment(segid='RNAA', pdbfile='chainB_cg.pdb')
gen.read_coords(segid='RNAA', filename='chainB_cg.pdb')

gen.write_pdb(filename='1jid.pdb')