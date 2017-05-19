#===========================================================================================
# This example read a density map and performs operations such as translation, rotation, etc. 
#===========================================================================================

from TEMPy.MapParser import MapParser
from TEMPy.StructureParser import PDBParser
import os



path_out='Test_Files'
if os.path.exists(path_out)==True:
	print "%s exists" %path_out
else:
	os.mkdir(path_out)
os.chdir(path_out)

#Example of user defined parameters
sim_res=6.6 #Target resolution of the outputted map.
sim_sigma_coeff=0.187 #Sigma width of the Gaussian used to blur the atomic structure.

structure_instance=PDBParser.read_PDB_file('1J6Z','1J6Z.pdb',hetatm=False,water=False)
print structure_instance

map_target=MapParser.readMRC('emd_5168_monomer.mrc') #read target map
print map_target
print map_target.get_com()
print "map information"
print "x"
print map_target.x_origin()
print "y"
print map_target.y_origin()
print "z"
print map_target.z_origin()
print "box size"
print map_target.box_size()	
print "map_size"
print map_target.map_size()
print 'min and max map'
print map_target.min()
print map_target.max()
print 'scale_map'
print map_target.scale_map(10)

print "Modifying maps"
print "rotate map"
map_target2 = map_target.rotate_by_axis_angle(10,3.7,12.4, 45, structure_instance.CoM) # Rotation around centre of mass of structure_instance
print map_target2
print "translate map"
map_target3 = map_target.translate(5.2,5,1) # Translation. Uses fourier-shifting, so movements are periodic.
print map_target3
print "normalise map"
map_target4 = map_target.normalise() # Normalising (mean=0, sd=1)
print map_target4
