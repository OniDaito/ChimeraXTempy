#===============================================================
# This example reads a PDB file and creates a structure instance 
#===============================================================

from TEMPy.StructureParser import PDBParser
import os


path_out='Test_Files'
if os.path.exists(path_out)==True:
	print "%s exists" %path_out
else:
	os.mkdir(path_out)
os.chdir(path_out)

#Generate Structure Instance from PDB File, default hetatm=False and water= False
'fetch a structure PDB file and create a structure instance'
structure_instance=PDBParser.fetch_PDB('1A5T','1A5T.pdb',hetatm=True,water=False)
print structure_instance

'fetch a structure mmCIF file and create a structure instance'
#need last version Biopython (1.40b)
#structure_instance=mmCIFParser.fetch_mmCIF('1A5T','1A5T.cif',hetatm=True,water=True)

'read a PDB files and create a structure instance' 
structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)
print structure_instance