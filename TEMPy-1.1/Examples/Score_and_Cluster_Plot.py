#==================================================================================================================    
# This example performs ranking by CCC and RMSD hierarchical clustering on an ensemble and plots it as a dendrogram
#==================================================================================================================    

from TEMPy.StructureParser import PDBParser  
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.EnsembleGeneration import  EnsembleGeneration
from TEMPy.Cluster import  Cluster
from TEMPy.ShowPlot import Plot
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
sim_map = blurrer.gaussian_blur(structure_instance, 6.6,densMap=map_target)

#list_rotate_models=EnsembleGeneration.randomise_structs(structure_instance, 20, 10, 60, v_grain=30, rad=False,write=False)
path_dir1="CCC_top20"
list_rotate_models=EnsembleGeneration.loadEnsemble(path_dir1,'mod_',hetatm=False,water=False,pdb=True)
list_rotate_models.append(['1J6Z',structure_instance])
Cluster=Cluster()

ranked_ensemble=Cluster.rank_fit_ensemble(list_rotate_models,'CCC',6.6,0.187,number_top_mod=20,write=False,targetMap=map_target.copy(),cont_targetMap=10.0)
mxRMSD=Cluster.RMSD_ensemble(ranked_ensemble,list_rotate_models)

#Plot 
cutoff =mxRMSD.mean()
Plot=Plot()
cluster_output=Plot.ShowHierarchicalClusterings(ranked_ensemble,mxRMSD,cutoff,name='CCCtop20_cluster',save=False,cluster_index=True)          
Plot.PrintOutClusterAnalysis(cluster_output)
