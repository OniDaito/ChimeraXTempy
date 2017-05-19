#===========================================================================================
# This example reads a structure and performs operations such as translation, rotation, etc. 
#===========================================================================================

from TEMPy.StructureParser import PDBParser
import os



path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)

#working dir    
os.chdir(path_out)


print 'read a PDB files and create a structure instance' 
structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=True,water=False)
print structure_instance

 
print "Translating a Structure (moved +4.3A in the x-direction, 1A in y and -55A in z)"
print "starting structure:"
print structure_instance
print "Centre of mass:"
print structure_instance.CoM 
print "translate the structure instance" 
structure_instance.translate(4.3, 1.0, -55)

structure_instance.write_to_PDB('testA.pdb')
print "Centre of mass of translated structure instance:"
print structure_instance.CoM
print "RMSD from initial position after translation"
print structure_instance.RMSD_from_init_position()
print "NOTE reset transformation with reset_position():"
structure_instance.reset_position() #note reset transformation.
  
print "get_vector_list"
print "Return an array containing Vector instances containing 3D coordinates of all the atoms"
print "here the position of the first atom:"
print structure_instance.get_vector_list()[0]

print "get_pos_mass_list"
print "Return an array containing Vector instances containing 3D coordinates of all atoms and mass"
print structure_instance.get_pos_mass_list()[0]
 
print "get_extreme_values"
print """Return a 6-ple containing the minimum and maximum x, y and z co-ordinates of the structure.
Given in order (min_x, max_x, min_y, max_y, min_z, max_z)."""
print structure_instance.get_extreme_values()
 
print "get_atom_list"
print "Return an array containing Atom instances of positions of all atoms."
print "here the position of the first atom:"
print structure_instance.get_atom_list()[0]
 
print "get_chain_list"
print "Return list of chain ID in Structure"
print structure_instance.get_chain_list()

print "get_chain"
print "return specific chain in structure"
print "structure instance chain A:"
print structure_instance.get_chain("A")

print "structure instance selection residues 67-160:"
print structure_instance.get_selection('7','160')
 
print 'break_into_segments used to create a Rigid Body'
rigid_list=[[151,156],[161,167]]# select two segments: 1st from residues 130 to 166, 2nd from residues 235 to 280
list_structure=structure_instance.break_into_segments(rigid_list)
print "here the list of the segments:"
print list_structure
 
print "add the selected segments to an existing Structure Instance"
print "the existing Structure Instance:"
print structure_instance
print "the modified one:"
print structure_instance.combine_structures(list_structure)
print "NOTE: the residues are NOT renumbered." 

print "combine a list of selected segments in a unique Structure Instance (rigid body)"
print "the starting Structure Instance:"
print structure_instance
print "selected rigid body:"
rigid_body1=structure_instance.combine_SSE_structures(list_structure)# create a Structure Instance of selected segments.
print rigid_body1

print 'get_selection_more_than'
print "from residues 167 till the end of the existing Structure Instance"
print structure_instance.get_selection_more_than('167') #selet from residues 167 till end
 
print 'get_residue'
print "return all atoms in a specified residues:"
print structure_instance.get_residue('167')
 
print 'get_atom'
print "return atom information"
print structure_instance.get_atom(200)
 
print 'get_backbone'
print "return structure instance that has only backbone atoms"
print structure_instance.get_backbone()
 
