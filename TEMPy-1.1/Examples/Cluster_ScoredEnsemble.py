#====================================================================
# This example performs hierarchical clustering based on CCC and RMSD 
#====================================================================

from TEMPy.StructureParser import PDBParser  
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.EnsembleGeneration import  EnsembleGeneration
from TEMPy.Cluster import  Cluster
import os

path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)
os.chdir(path_out)

structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)

blurrer = StructureBlurrer()
EnsembleGeneration=EnsembleGeneration()
scorer = ScoringFunctions()

map_target=MapParser.readMRC('emd_5168_monomer.mrc') #read target map
print map_target

map_probe = blurrer.gaussian_blur(structure_instance, 6.6,densMap=map_target)
list_rotate_models=EnsembleGeneration.randomise_structs(structure_instance, 20, 10, 60, v_grain=30, rad=False,write=False)

Cluster=Cluster()
ranked_ensemble=Cluster.cluster_fit_ensemble_top_fit(list_rotate_models,'CCC',1.5,6.6,0.187,number_top_mod=3,write=False,targetMap=map_target.copy())
 