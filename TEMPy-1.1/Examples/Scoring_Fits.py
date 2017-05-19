#============================================================
# This example scores a fit using different scoring functions 
#============================================================

from TEMPy.MapParser import MapParser
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
import os


path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)
os.chdir(path_out)

#read PDB file and create a Structure instance.
#note hetatm and water to include
structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)
print "structure_instance:"
print structure_instance

blurrer = StructureBlurrer()
scorer = ScoringFunctions()

map_target=MapParser.readMRC('emd_5168_monomer.mrc') #read target map

map_probe = blurrer.gaussian_blur(structure_instance, 6.6,densMap=map_target)#create a simulated map from the structure instance
map_probe.write_to_MRC_file("map_probe_actin.mrc") #write simulated map to a MRC file format

##SCORING FUNCTION

print "Calculate Envelope Score (ENV):"
molecualr_weight=structure_instance.get_prot_mass_from_atoms()
#Mmolecualr_weight=structure_instance.get_prot_mass_from_res()
first_bound=map_target.get_primary_boundary(molecualr_weight, map_target.min(), map_target.max())
#print scorer.envelope_score_APJ(map_target, first_bound, structure_instance,norm=True)
print scorer.envelope_score(map_target, first_bound, structure_instance,norm=True)

print "Calculate Mutual information Score (MI)"
print scorer.MI(map_target,map_probe)

print "Calculate Laplacian cross-correlation Score (LAP)"
print scorer.laplace_CCC(map_target,map_probe)
 
print "Calculate cross-correlation Score (CCC)"
print scorer.CCC(map_target,map_probe)

print "Calculate Normal Vector (NV):"
#Number of points to use in the normal vector score
points= round((map_target.map_size())*0.01)
first_bound=map_target.get_primary_boundary(structure_instance.get_prot_mass_from_atoms(), map_target.min(), map_target.max())
second_bound=map_target.get_second_boundary(first_bound, points, first_bound, map_target.max(),err_percent=1)
print scorer.normal_vector_score(map_target,map_probe, first_bound, second_bound)

print "Calculate Normal Vector (NV) with Sobel Filter:"
print scorer.normal_vector_score(map_target,map_probe, first_bound, second_bound,Filter='Sobel')
