#========================================================================
# This example generates a random ensemble of fits and scores it with CCC 
#========================================================================

from TEMPy.StructureParser import PDBParser  
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.EnsembleGeneration import  EnsembleGeneration
import os

path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)
os.chdir(path_out)


structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)
print structure_instance

blurrer = StructureBlurrer()
EnsembleGeneration=EnsembleGeneration()
scorer = ScoringFunctions()

map_target=MapParser.readMRC('emd_5168_monomer.mrc') #read target map
map_probe = blurrer.gaussian_blur(structure_instance, 6.6,densMap=map_target)#create a simulated map from the structure instance

#Create a Random ensemble of 10 structures randomly within  5 A translation and 60 deg rotation.
list_rotate_models=EnsembleGeneration.randomise_structs(structure_instance, 10, 5, 60, v_grain=30, rad=False,write=True)


#CCC score from starting fit
line='%s %s\n'%('1J6Z',scorer.CCC(map_probe,map_target))
count=0
#loop to score each of the alternative fits in the ensemble
for mod in list_rotate_models:
        count+=1
        mod_name=mod[0]
        mod_structure_instance=mod[1]
        map_probe = blurrer.gaussian_blur(mod_structure_instance, 6.6,densMap=map_target,sigma_coeff=0.187)
        line+='%s %s\n'%(mod_name,scorer.CCC(map_probe,map_target))        
print line
        
