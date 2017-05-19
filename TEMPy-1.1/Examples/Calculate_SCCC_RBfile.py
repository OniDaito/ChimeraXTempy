#=============================================================================================
# This example calcualtes SCCC for local fit on a given set of rigid bodies (read from a file)
#=============================================================================================

from TEMPy.StructureParser import PDBParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.ShowPlot import Plot
import os

sim_sigma_coeff=0.187 #Sigma width of the Gaussian used to blur the atomic structure.

path_out='Test_Files/'
if os.path.exists(path_out)==True:
	print "%s exists" %path_out
else:
	os.mkdir(path_out)
os.chdir(path_out)

rb_file ="3MFP_RigidBodyFile.txt" #Rigid-body file alike in which the selection of each segments is specified. To use for the SSCCC score.


#rb_file2 ="1J6Z_sse.txt"

structure_instance=PDBParser.read_PDB_file('3MFP','3MFP.pdb',hetatm=False,water=False)
print structure_instance

structure_instance2=PDBParser.read_PDB_file('1J6Z.pdb','1J6Z.pdb',hetatm=False,water=False)
print structure_instance2


blurrer = StructureBlurrer()
scorer = ScoringFunctions()
Plot=Plot()

emmap=MapParser.readMRC('emd_5168_monomer.mrc') #read target map
print emmap

sim_map = blurrer.gaussian_blur(structure_instance, 6.6,densMap=emmap,sigma_coeff=sim_sigma_coeff,normalise=True)
print 'structure_instance',scorer.CCC(sim_map,emmap)
print sim_map


sim_map2 = blurrer.gaussian_blur(structure_instance2, 6.6,densMap=emmap,sigma_coeff=sim_sigma_coeff,normalise=True)
print 'structure_instance_same',scorer.CCC(sim_map2,emmap)

SCCC_list_structure_instance=[]
listRB=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance2)
for RB in listRB:
		score_SCCC=scorer.SCCC(emmap,6.6,sim_sigma_coeff,structure_instance2,RB)
		SCCC_list_structure_instance.append(score_SCCC)
		print score_SCCC

listRB=RBParser.RBfileToRBlist(rb_file)
Plot.PrintOutChimeraAttributeFileSCCC_Score('3MFP',SCCC_list_structure_instance,listRB)


SCCC_list_structure_instance2=[]
listRB2=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance2)

for RB in listRB2:
		score_SCCC=scorer.SCCC(emmap,6.6,sim_sigma_coeff,structure_instance2,RB)
		SCCC_list_structure_instance2.append(score_SCCC)
		
listRB2=RBParser.RBfileToRBlist(rb_file)
Plot.PrintOutChimeraAttributeFileSCCC_Score('1J6Z',SCCC_list_structure_instance2,listRB2)


