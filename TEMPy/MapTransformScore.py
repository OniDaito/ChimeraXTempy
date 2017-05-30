#===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#	  Copyright  2015 Birkbeck College University of London. 
#
#				Authors: Maya Topf, Daven Vasishtan, Arun Prasad Pandurangan,
#						Irene Farabella, Agnel-Praveen Joseph, Harpal Sahota
# 
#     This software is made available under GPL V3 license
#     http://www.gnu.org/licenses/gpl-3.0.html
#     
#     
#     Please cite your use of TEMPy in published work:
#     
#     Farabella, I., Vasishtan, D., Joseph, A.P., Pandurangan, A.P., Sahota, H. & Topf, M. (2015). J. Appl. Cryst. 48.
#
#===============================================================================

from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions

class MapTransformationScore:
    """A class to create ensemble and score of Map instance"""
    
    def __init__(self):
        pass
       
    def transform_map(self,matR, transvec, m1, m2, c1, c2):
        mat = matR.T
        emmap1=MapParser.readMRC(m1)
        emmap2=MapParser.readMRC(m2)
        # geometric centre of map
        vec_centre = emmap2.centre()
        spacing = emmap2.apix
        # to work on the box transformations, get the box centre irrespective of origin
        vec_centre.x = vec_centre.x - emmap2.x_origin()
        vec_centre.y = vec_centre.y - emmap2.y_origin()
        vec_centre.z = vec_centre.z - emmap2.z_origin()

        # calculate new box dimensions, after rotation
        new_centre = emmap2._box_transform(matR)

        output_shape = (int(new_centre.x/spacing),int(new_centre.y/spacing),int(new_centre.z/spacing))
        new_centre.x = new_centre.x/2
        new_centre.y = new_centre.y/2
        new_centre.z = new_centre.z/2
        # offset for rotation
        offset = emmap2._rotation_offset(mat, vec_centre, new_centre)

        #APPLY ROTATION
        emmap2 = emmap2._matrix_transform_offset(mat, output_shape, offset.x, offset.y, offset.z)

        offset_x = new_centre.x - vec_centre.x
        offset_y = new_centre.y - vec_centre.y
        offset_z = new_centre.z - vec_centre.z
        emmap2 = emmap2.shift_origin(-offset_x,-offset_y,-offset_z)

        # TRANSLATION COMPONENT
        a14,a24,a34 =  transvec[0],transvec[1],transvec[2]
        emmap_2 = emmap2.shift_origin(float(a14)*spacing,float(a24)*spacing,float(a34)*spacing)

        emmap_1 = emmap1.copy()
        # CROP BOX TO REDUCE ARRAY SIZE
        emmap_1._crop_box(c1,2)
        emmap_2._crop_box(c2,2)


        # DETERMINE A COMMON ALIGNMENT BOX
        spacing = emmap_2.apix
        if emmap_2.apix < emmap_1.apix: spacing = emmap_1.apix
        grid_shape, new_ori = emmap_1._alignment_box(emmap_2,spacing)

        # INTERPOLATE TO NEW GRID
        emmap_1 = emmap_1._interpolate_to_grid(grid_shape,spacing,new_ori)
        emmap_2 = emmap_2._interpolate_to_grid(grid_shape,spacing,new_ori)

        sc = ScoringFunctions()
        ccc = sc.CCF_mask_zero(emmap_1,emmap_2,c1,c2)
        mi = sc.MI(emmap_1,emmap_2)
        env = sc.map_envelope_score(emmap_1,emmap_2,c1,c2)
        nv = sc.normal_vector_score(emmap_1,emmap_2,float(c1)-(emmap1.std()*0.05),float(c1)+(emmap1.std()*0.05))
        nv = sc.normal_vector_score(emmap_1,emmap_2,float(c1)-(emmap1.std()*0.05),float(c1)+(emmap1.std()*0.05),Filter='Sobel')

        return ccc, mi, env, nv ,nv_s
