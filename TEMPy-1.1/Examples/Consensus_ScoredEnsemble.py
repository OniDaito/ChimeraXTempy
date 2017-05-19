#==================================================================================================================    
# This example performs ranking by CCC and RMSD hierarchical clustering on an ensemble and plots it as a dendrogram
#==================================================================================================================    

from TEMPy.StructureParser import PDBParser  
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.EnsembleGeneration import  EnsembleGeneration
from TEMPy.Consensus import Consensus
import os



path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)

os.chdir(path_out)

structure_instance=PDBParser.read_PDB_file('1GGG','1ggg.pdb',hetatm=False,water=False)

blurrer = StructureBlurrer()
EnsembleGeneration=EnsembleGeneration()
scorer = ScoringFunctions()
map_target=MapParser.readMRC('1ggg_5A.mrc') #read target map
sim_map = blurrer.gaussian_blur(structure_instance, 6.6,densMap=map_target)

#list_rotate_models=EnsembleGeneration.randomise_structs(structure_instance, 20, 10, 60, v_grain=30, rad=False,write=False)
path_dir1="CCC_top20"
list_rotate_models=EnsembleGeneration.randomise_structs(structure_instance, 5, 10, 60, v_grain=30, rad=False,write=False)
Consensus=Consensus()
score_list=['CCC','MI','NV_Sobel','LAP','ENV']
Consensus.vote(list_rotate_models,score_list,6.6,0.187,number_top_mod=0,targetMap=map_target.copy())
