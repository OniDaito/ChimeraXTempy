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

from numpy import linspace
import os

import TEMPy.Vector as Vector
from TEMPy.StructureParser import PDBParser


class EnsembleGeneration:
    """A class to create ensemble of structure instance"""
    
    def __init__(self):
        pass
    
    def loadEnsemble(self,path_dir,file_name_flag,hetatm=False,water=False,verbose=False,pdb=True):
        """
        Load an ensemble of Structure Instance from the directory path_dir.
        Arguments:
            *path_dir*
                directory name
            *file_name_flag*
                name or suffix of the files.
        """
        structure_list=[]
        list_rotate_models=[filein for filein in os.listdir(path_dir) if file_name_flag in filein and filein[-4:]=='.pdb' ]
        for pdbin in list_rotate_models:
            print pdbin
            if pdb==True:
                file_in=path_dir+'/'+pdbin
                #print file_in
                if verbose==True:
                    print "load file:",pdbin[:-4],file_in
                structure_instance=PDBParser.read_PDB_file(str(pdbin[:-4]),str(file_in),hetatm=hetatm,water=water)
                structure_list.append([pdbin[:-4],structure_instance])
        return structure_list

    def randomise_structs(self,structure_instance, no_of_structs, max_trans, max_rot, v_grain=30, rad=False,flag='mod',write=False):
        """
          
        Generate an ensemble of Structure Instance. 
                     
        Arguments:
            *structure_instance*
                Input Structure Instance
            *no_of_structs*
                int, number of structures to output
            *max_trans*
                Maximum translation permitted
            *max_rot*
                Maximum rotation permitted (in degree if rad=False)
            *v_grain*
                Graning Level for the generation of random vectors (default=30) 
            *write*
                True will write out each structure_instanceure Instance in the ensemble as single PDB.
        Return:
            list of Structure Instance in which each item is [structure_instance_name,structure_instance]
        """
        ensemble_list=[]
        file_name0=flag+'_0'
        ensemble_list.append([file_name0,structure_instance.copy()])
        if write==True:
            structure_instance.write_to_PDB(file_name0)
        count=0
        for x in range(0,no_of_structs-1):
            count+=1
            file_name=flag+'_'+str(count)+".pdb"
            structure_instance.randomise_position(float(max_trans), float(max_rot), v_grain, rad)
            if write==True:
                structure_instance.write_to_PDB(file_name)
            ensemble_list.append([file_name[:-4],structure_instance.copy()])
            structure_instance.reset_position() 
        return ensemble_list
   
    def anglar_sweep(self,structure_instance, axis, translation_vector, no_of_structs, loc_rotation_angle, flag='mod', atom_com_ind=False,write=False,filename=None):
        """
        Generate an ensemble of Structure Instance
        
        NOTE - Chose the number of structures for the ensemble accordingly with the angular increment step (loc_rotation_angle/no_of_structs) and
        translational increment step (translation_vector/no_of_structs) required. 
        Default setting is around the center of mass.
        
        Arguments:
            *structure_instance*
                Input Structure Instructure_instance
            *axis*
                3-tuple, axis for translation
            *translation_vector*
                3-ple, vector for translation
            *no_of_structs*
                int, number of structures to output
            *loc_rotation_angle*
                tuple, rotation angle for local rotation (degrees)
            *flag*
                string, prefix name for outputted pdb files
            *atom_com_ind*
                int, index of atom to rotate around. If False, rotates around centre of mass
            *write*
                True will write out each Structure Instance in the ensemble as single PDB.
        Return:
            list of Structure Instance in which each item is [structure_instance_name,structure_instance]
        """
        
            # Work out distance between adjacent structures
        ensemble_list=[]
        
        file_name0=flag+'_0'
        ensemble_list.append([file_name0,structure_instance])
        
        grain = loc_rotation_angle/no_of_structs
        if int(grain)<1:
            print "Warning: less then 1deg rotation"
        else:
            transl_x=linspace(0, translation_vector[0], num=no_of_structs)
            transl_y=linspace(0, translation_vector[1], num=no_of_structs)
            transl_z=linspace(0, translation_vector[2], num=no_of_structs)
            angle_rot=linspace(0,loc_rotation_angle,num=no_of_structs)
            count=-1    
            for x in angle_rot:
                count+=1
                file_name=flag+'_'+str(count+1)+".pdb"
                if atom_com_ind:
                    loc_point = structure_instance[atom_com_ind].get_pos_vector()
                    structure_instance.rotate_by_axis_angle(axis[0],axis[1],axis[2], x, com=loc_point)
                else:
                    structure_instance.rotate_by_axis_angle(axis[0],axis[1],axis[2],x)
                structure_instance.translate(transl_x[count],transl_y[count],transl_z[count])
                if write==True:
                    #print "model_num: ",(count+1),"angle: ",x,"tranls_increment: ",transl_x[count],transl_y[count],transl_z[count]
                    structure_instance.write_to_PDB(file_name)
                    ensemble_list.append([file_name[:-4],structure_instance])
                else:
                    ensemble_list.append([file_name[:-4],structure_instance])
                    #print "model_num: ",(count+1),"angle: ",x,"tranls_increment: ",transl_x[count],transl_y[count],transl_z[count]
                    structure_instance.reset_position()
        return ensemble_list
 
    def spiral_sweep(self, structure_instance, axis, dist, no_of_structs, loc_ang_range, loc_axis, flag, atom_com_ind=False,write=False):
        """
        Generate an ensemble of Structure Instance
        
        Arguments:
            *structure_instance*
                Input Structure Instance
           *axis*
               3-ple, axis for translation
           *dist*
               int, translation range (Angstroms)
           *no_of_structs*
               int, number of structures to output
           *loc_axis*
               3-ple, axis for local rotation around centre of mass
           *loc_ang_range*
                tuple, rotation range for local rotation (degrees)
           *flag*
               string, prefix name for outputted pdb files
           *atom_com_ind*
               int, index of atom to rotate around. If False, rotates around centre of mass.
            *write*
                True will write out each Structure Instance in the ensemble as single PDB.
         Return:
            list of Structure Instance in which each item is [structure_instance_name,structure_instance]
               
        """
        ensemble_list=[]
        file_name0='mod_0'
        ensemble_list.append([file_name0,structure_instance])
        count=0
        # Work out distance between adjacent structures
        grain = dist/no_of_structs
        # Make axis into vector of length 1
        axis = Vector.Vector(axis[0], axis[1], axis[2]).unit()

        for r in range(no_of_structs):
            file_name=flag+str(r+1)+".pdb"
            # Translate structure along axis
            structure_instance.translate(axis.x*r*grain,axis.y*r*grain,axis.z*r*grain)

            # Rotation around centre of mass
            loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
            if atom_com_ind:
                loc_point = self[atom_com_ind].get_pos_vector()
                structure_instance.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain, com=loc_point)
            else:
                structure_instance.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)

                if write==True:

                    structure_instance.write_to_PDB(file_name)
                    ensemble_list.append([file_name[:-4],structure_instance])
                else:
                    ensemble_list.append([file_name[:-4],structure_instance])
                structure_instance.reset_position()
        return ensemble_list


#TO ADD LOAD ENSEMBLE
#===============================================================================
# def load_dir_pdb(path_dir,file_name_flag,hetatm=False,water=False):
#     
#     structure_list=[]
#     list_rotate_models=[file for file in os.listdir(path_dir) if file_name_flag in file]
#     for pdb in list_rotate_models:
#         file_in=path_dir+'/'+pdb
#         print "*****",pdb[:-4],file_in
#         structure_instance=PDBParser.read_PDB_file(str(pdb[:-4]),str(file_in),hetatm=hetatm,water=water)
#         structure_list.append([pdb[:-4],structure_instance])
#     return structure_list
#===============================================================================
                
#NEED TESTING
#
#
#===============================================================================
#     def circular_sweep_with_fixed_spin(self,struct, axis, point, ang_range, no_of_structs, flag, loc_axis=True, loc_ang_range=[45,90], atom_com_ind=False,write=False):
#         """
#         Generate an ensemble of Structure Instance
#         
#         Arguments:
#             *struct*
#                 Input Structure Instance
#             *axis*
#                 3-ple, axis around which the large scale rotation will be done
#             *point*
#                 3-ple, point around which large scale rotation will be done
#             *ang_range*
#                 tuple, rotation range (degrees)
#             *no_of_structs*
#                 int, number of structures to output
#             *flag*
#                 string, prefix name for outputted pdb files
#             *loc_axis*
#                 3-ple, axis for local rotation around centre of mass
#            *loc_ang_range*
#                tuple, rotation range for local rotation (degrees)
#            *atom_com_ind*
#                int, index of atom to rotate around. If False, rotates around centre of mass.
#             *write*
#                True will write out each Structure Instance in the ensemble as single PDB.
#                 
#         """
#         ensemble_list=[]
#         file_name0='mod_0'
#         ensemble_list.append([file_name0,struct])
#         # Work out angle between adjacent structures
#         grain = (ang_range[1]-ang_range[0])/no_of_structs
#         # Make point into Vector object
#         point = Vector(point[0], point[1],point[2])
# 
#         for r in range(no_of_structs):
#             # Rotation around defined point
#             struct.rotate_by_axis_angle(axis[0],axis[1],axis[2], r*grain, com=point)
#             print r
#             if loc_axis and loc_ang_range:
#                 # Rotation around centre of mass or atom, if specified
#                 loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
#                 if atom_com_ind:
#                     loc_point = a[atom_com_ind].get_pos_vector()
#                     struct.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain, com=loc_point)
#                 else:
#                     struct.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)
#                 if write==True:
#                     struct.write_to_PDB(file_name)
#                     ensemble_list.append([file_name[:-4],struct])
#                 else:
#                     ensemble_list.append([file_name[:-4],struct])
#                     print "model_num: ",(count+1),"angle: ",x,"tranls_increment: ",transl_x[count],transl_y[count],transl_z[count]
#                     struct.reset_position()
#             else:
#                 print "print joke"
#         return ensemble_list
# 
# 
# 
# 
#     def circular_sweep_with_tangential_spin(self,struct, axis, point, ang_range, no_of_structs, flag, loc_ang_range=False, atom_com_ind=False,write=False):
#  
#         """
#         Generate an ensemble of Structure Instance
#         
#         Arguments:
#             *struct*
#                 Input Structure Instance
#             *axis*
#                 3-ple, axis around which the large scale rotation will be done
#             *point*
#                 3-ple, point around which large scale rotation will be done
#             *ang_range*
#                 tuple, rotation range (degrees)
#             *no_of_structs*
#                 int, number of structures to output
#             *flag*
#                 string, prefix name for outputted pdb files
#             *loc_ang_range*
#                tuple, rotation range for local rotation (degrees)
#             *atom_com_ind*
#                 int, index of atom to rotate around. If False, rotates around centre of mass.
#             *write*
#                 True will write out each Structure Instance in the ensemble as single PDB.
#                 
#         """
#         ensemble_list=[]
#         file_name0='mod_0'
#         ensemble_list.append([file_name0,struct])
#         
#         # Work out angle between adjacent structures
#         grain = (ang_range[1]-ang_range[0])/no_of_structs
#         # Make point and axis into Vector objects
#         point = Vector(point[0], point[1],point[2])
#         axis = Vector(axis[0], axis[1], axis[2])
# 
#         for r in range(no_of_structs):
#             # Rotation around defined point
#             self.rotate_by_axis_angle(axis[0],axis[1],axis[2], r*grain, com=point)
# 
#             if loc_ang_range:
#                 # Rotation around centre of mass or atom, if specified
#                 loc_grain = (loc_ang_range[1]-loc_ang_range[0])/no_of_structs
#                 if atom_com_ind:
#                     loc_point = self[atom_com_ind].get_pos_vector()
#                 else:
#                     loc_point = self.CoM    
#                 # Axis of local spin is cross product of axis and radius of rotation (ie. the tangent of the circular arc being traversed)
#                 rad = (loc_point-point).unit()
#                 loc_axis = rad.cross(axis)
#                 self.rotate_by_axis_angle(loc_axis[0],loc_axis[1],loc_axis[2], r*loc_grain)
# 
#                 if write==True:
#                     #print "model_num: ",(count+1),"angle: ",x,"tranls_increment: ",transl_x[count],transl_y[count],transl_z[count]
#                     struct.write_to_PDB(file_name)
#                     ensemble_list.append([file_name[:-4],struct])
#                 else:
#                     ensemble_list.append([file_name[:-4],struct])
#                     #print "model_num: ",(count+1),"angle: ",x,"tranls_increment: ",transl_x[count],transl_y[count],transl_z[count]
#                     struct.reset_position()
#         return ensemble_list
# 
#         
#===============================================================================

