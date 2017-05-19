#========================================================================================================================
# This example performs exhaustive map trasformation from an initial aligned pair of maps and calcualtes CCC, MI, ENV, NV
#========================================================================================================================

import os
import sys
from TEMPy.Vector import euler_to_matrix
from TEMPy.MapTransformScore import MapTransformationScore
import random

path_out='Test_Files'
if os.path.exists(path_out)==True:
    print "%s exists" %path_out
else:
    os.mkdir(path_out)
os.chdir(path_out)

# Rotate by 30 degree steps 
# Note: this example does an exhaustive search. We recommend to use a 10 degree rotation step for a more detailed sampling.
rot = 30
list_angles = range(360/rot)
# translate by 1 voxel steps
trans = 1
list_moves = range(20)
m1 = os.getcwd()+'/emd_5170.map'
m2 = os.getcwd()+'/emd_5168.map'
#recommend contour level for the map by the authors 
c1 = 0.61
c2 = 10.0

list_matrices = []
list_ccc = []
list_mi = []
list_env = []
list_nv=[]
list_nvs = []
list_transvec = []

index_trans = 0
MapTransformationScore=MapTransformationScore()

for i in random.sample(list_angles,5):
    mat1 = euler_to_matrix(i*rot,0,0)
    matR = mat1
    for j in random.sample(list_moves,3):
        listrans = [0.0,0.0,0.0]
        for p in random.sample(range(3),1):
            listrans[p] = j*trans
            ccc, mi, env, nv, nvs = MapTransformationScore.transform_map(matR, listrans, m1, m2, c1, c2)
            list_ccc.append(ccc)
            list_mi.append(mi)
            list_env.append(env)
            list_nv.append(nv)
            list_nvs.append(nvs)
            index_trans += 1
            list_transvec.append(listrans)
            list_matrices.append(mat1)
            print 'Rotation (degrees): ', i*rot
            print "Scores for translation %s is ccc: %f,mi: %f,env: %f,nv: %f,nvs: %f"%(listrans,ccc,mi,env,nv,nvs)
	    if j == 0: break		
            listrans = [0.0,0.0,0.0]

listcc = list_ccc[:]
list_ccc.sort()
ind = listcc.index(list_ccc[-1])
print '### Transformation for best local cross correlation score ###'
print 'Rotation Matrix is : ', list_matrices[ind]
print 'Translation Vector : ', list_transvec[ind]
print 'Cross correlation score : ', listcc[ind]
print 'Mutual Information score : ', list_mi[ind]
print 'Envelope score : ', list_env[ind]
print 'Normal Vector Score : ', list_nv[ind]
print 'Normal Vector Score with Sobel Filter : ', list_nvs[ind]

