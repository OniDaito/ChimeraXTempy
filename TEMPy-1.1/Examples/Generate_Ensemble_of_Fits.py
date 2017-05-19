#===============================================================================================
# This example generates an ensemble of configurations from the initial fit using angular sweeps
#===============================================================================================

from TEMPy.StructureParser import PDBParser
from TEMPy.EnsembleGeneration import  EnsembleGeneration
import os



path_out='Test_Files'
if os.path.exists(path_out)==True:
	print "%s exists" %path_out
else:
	os.mkdir(path_out)
os.chdir(path_out)


'read a PDB files and create a structure instance' 
structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)
print structure_instance


translation_vector=[4.3, 1.0, -55]
rotation_angle= 110
axis=[0.21949010788898163, -0.80559787935161753, -0.55030527207975843]
print "rotation: ",rotation_angle
print "axix: ",axis
print "translation_vector",translation_vector
print "generate angular sweep for 1J6Z"
#EnsembleGeneration=EnsembleGeneration()
list_ensemble=EnsembleGeneration().anglar_sweep(structure_instance,axis, translation_vector, 10, rotation_angle, 'mdl_angular_sweep', atom_com_ind=False)


for struct in list_ensemble:
	print "structure name: ",struct[0],struct[1]
