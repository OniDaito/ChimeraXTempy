#===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#	  Copyright 2015 Birkbeck College University of London. 
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

from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.ScoringFunctions import ScoringFunctions
from numpy import zeros
import sys

class Cluster:
    """A class to clustering an ensemble of structure instance"""
    
    def __init__(self):
        pass

    def _print_results_cluster(self,models,class_num,number_top_mod,score,write=False):
        """
        private function used in Cluster_Ensemble
        """
        out_list=[]
        if write==True:
            outp = open("top"+str(number_top_mod)+str(score)+"_classes.txt", "w")
            outp.write("pdb_name\tscore\tlrms\tclass\n")

            for i in range(1,class_num+1):
       
    # print the fits of each class ordered by the highest score 
              for ipdb in models:
                 if (ipdb[-1] == i):
                    out_list.append([ipdb[0],ipdb[2],ipdb[3],ipdb[4]])
                    outp.write("%s\t%.5f\t%.3f\t%d\n" %(ipdb[0],ipdb[2],ipdb[3],ipdb[4]))
            outp.close()
        else:
            for i in range(1,class_num+1):
              for ipdb in models:
                 if (ipdb[-1] == i):
                    out_list.append([ipdb[0],ipdb[2],ipdb[3],ipdb[4]])
        return out_list

    def _print_results_cluster2(self,models,write=True):
        """
        private function used in Cluster_Ensemble
        """
        out_list=[]
        if write==True:
            outp = open("top_rank.txt", "w")
            outp.write("pdb_name\tscore\tlrms\tclass\n")

            for i in models:
               #[name_mod,mod,score_mod,int(0),int(0)]
    # print the fits of each class ordered by the highest score 
                outp.write("%s\t%.5f\n" %(i[0],i[2]))
            outp.close()
        else:
            print 'this is for print!!!'

    def cluster_fit_ensemble_top_fit(self,ensemble_list,score,rms_cutoff,res_target_map,sigma_coeff,number_top_mod=0,write=False,targetMap=False):
        """
          
            RMSD clustering of the multiple "fits" starting from the best scoring model accordingly with a chosen score.
            Cluster the fits based on Calpha RMSD (starting from the best scoring model)            
                
                Arguments:
                    *ensemble_list*
                        Input list of Structure Instances.
                    *targetMap*
                        Target Map Instance.
                    *score*
                        Scoring function to use. 
                        See ScoringFunctions class for a list of the available Scoring Function.
                        E.g. set score='CCC' to use the Cross-correlation coefficient.
                        
                        Score option are:
                        
                        i    'CCC' - Cross-correlation coefficient; 
                        
                        ii    'LAP' - Laplacian-filtered cross-correlation coefficient:  useful for maps with resolutions worse than 10-15 A;
                        
                        iii   'MI' - Mutual information score: a good and robust score but relatively slow to calculate; 
                        
                        iv    'ENV' - Envelope score: the fastest score to calculate due to binarisation of the map. 
                        
                        v-vii  'NV','NV_Sobel','NV_Laplace'- Normal vector score: a vector-based surface superimposition score with or without Sobel/Laplace filter.

                        viii 'CD' - Chamfer Distance: a score used in computer vision algorithms as a fast similarity metric 

                    *rms_cutoff*
                        float,  the Calpha RMSD cutoff based on which you want to cluster the solutions. For example 3.5 (for 3.5 A).
                    *res_target_map*
                        the resolution, in Angstroms, of the target Map.
                    *sigma_coeff*
                        the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                        Default values is 0.356.
                
                        Other values used :
                
                            0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                            0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                            0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                            0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                            0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                            1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

                    *number_top_mod*
                        Number of Fits to cluster. Default is all.
                    *write*
                        True will write out a file that contains the list of the structure instances representing different fits scored and clustered.
                        note the lrms column is the Calpha RMSD of each fit from the first fit in its class
        """
        blurrer = StructureBlurrer()
        
        scorer = ScoringFunctions()
        
        cluster=Cluster()

        count=0
        dict_ensembl={}
        list_ordered=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=0,write=False,targetMap=targetMap.copy())
        #cluster fits by local rmsd
        if number_top_mod==0:
            ini_num = 0
            end_num = len(list_ordered)
            fit_class = 0
            for ipdb in list_ordered:
                print "model num %d: %s\n" %(list_ordered.index(ipdb)+1, ipdb[0])
                ini_num1 = list_ordered.index(ipdb)
                mod1=ipdb[1]
                print 'next index ' + str(ini_num1)
                if ipdb[-1] == 0:
                    fit_class+=1
                    for ipdb1 in list_ordered[ini_num1 : end_num]:
                        mod2=ipdb1[1]
                        if ipdb1[-1] == 0:
                            rmsd_val=float(mod1.RMSD_from_same_structure(mod2,CA=True))
                            ipdb1[3]=rmsd_val
                            print "rmsd of %s from best local fit (%s)= %.2f" %(ipdb1[0], ipdb[0], rmsd_val)
                            if rmsd_val < rms_cutoff:
                                ipdb1[-1] = fit_class
                            print 'class= ' + str(ipdb1[-1])
                        else: continue
                else: continue
            return cluster._print_results_cluster(list_ordered,fit_class,number_top_mod,score,write)
        else:
            x=int(number_top_mod)
            ini_num = 0
            end_num = len(list_ordered[:x])
            fit_class = 0
            for ipdb in list_ordered[:x]:
                print "model num %d: %s\n" %(list_ordered.index(ipdb)+1, ipdb[0])
                ini_num1 = list_ordered.index(ipdb)
                mod1=ipdb[1]
                print 'next index ' + str(ini_num1)
                if ipdb[-1] == 0:
                    fit_class+=1
                    for ipdb1 in list_ordered[ini_num1 : end_num]:
                        mod2=ipdb1[1]
                        if ipdb1[-1] == 0:
                            rmsd_val=float(mod1.RMSD_from_same_structure(mod2,CA=True))
                            print "rms of %s from best local fit (%s)= %.2f" %(ipdb1[0], ipdb[0], rmsd_val)
                            ipdb1[3]=rmsd_val
                            if rmsd_val < rms_cutoff:
                                ipdb1[-1] = fit_class
                            print 'class= ' + str(ipdb1[-1])
                        else: continue
                else: continue
            return cluster._print_results_cluster(list_ordered[:x],fit_class,number_top_mod,score,write)
        
    def RMSD_ensemble(self,rank_fit_ensemble,ensemble_list,CA=True):
        """
            Calculates the pairwise RMSD matrix for all Structure Instance in the ensemble.
        
                Arguments:
                    *rank_fit_ensemble*
                        Ensemble of Structure Instance ranked using cluster.rank_fit_ensemble
                    *ensemble_list*
                        Input list of Structure Instances
                    *CA is set to True if only CA-RMSD is needed*
                Return:
                    A numpy array
        """
        list_rotate_models_dict={}
        for i in ensemble_list:
            list_rotate_models_dict[i[0]]=i[1]
                    
        sorted_rank=rank_fit_ensemble
        mxRMSD = zeros(shape=(len(sorted_rank),len(sorted_rank)))
        for mod1 in sorted_rank:
            for mod2 in sorted_rank:
                print mod1[0],mod2[0]
                rmsd_val=float(list_rotate_models_dict[mod1[0]].RMSD_from_same_structure(list_rotate_models_dict[mod2[0]],CA=CA))
                m1=sorted_rank.index(mod1)
                m2=sorted_rank.index(mod2)
                mxRMSD[m1][m2]=rmsd_val
        return mxRMSD

    def rank_fit_ensemble(self,ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=0,\
                          write=False,targetMap=False,cont_targetMap=None):
        """
          
            RMSD clustering of the multiple "fits" accordingly with a chosen score.
            Cluster the fits based on Calpha RMSD (starting from the best scoring model)            
                Arguments:
                    *ensemble_list*
                        Input list of Structure Instances.
                    *targetMap*
                        Target Map Instance.
                    *score*
                        Scoring function to use. 
                        See ScoringFunctions class for a list of the available Scoring Function.
                        E.g. set score='CCC' to use the Cross-correlation coefficient.
                        
                        Score option are:
                        
                        i    'CCC' - Cross-correlation coefficient; 
                        
                        ii    'LAP' - Laplacian-filtered cross-correlation coefficient:  useful for maps with resolutions worse than 10-15 A;
                        
                        iii   'MI' - Mutual information score: a good and robust score but relatively slow to calculate; 
                        
                        iv    'ENV' - Envelope score: the fastest score to calculate due to binarisation of the map. 
                        
                        v-vii  'NV','NV_Sobel','NV_Laplace'- Normal vector score: a vector-based surface superimposition score with or without Sobel/Laplace filter.

                        viii 'CD' - Chamfer Distance: a score used in computer vision algorithms as a fast similarity metric 
                                                                                         

                    *rms_cutoff*
                        float,  the Calpha RMSD cutoff based on which you want to cluster the solutions. For example 3.5 (for 3.5 A).
                    *res_target_map*
                        the resolution, in Angstroms, of the target Map.
                    *sigma_coeff*
                        the sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                        Default values is 0.356.
                
                        Other values used :
                
                            0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                            0.225R which makes the Fourier transform of the distribution fall to 1/e of its maximum value at wavenumber 1/resolution, the default in Chimera (Petterson et al, 2004)
                    
                            0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, an option in Chimera (Petterson et al, 2004);
                    
                            0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                            0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                            1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).

                    *number_top_mod*
                        Number of Fits to cluster. Default is all.
                    *write*
                        True will write out a file that contains the list of the structure instances representing different fits scored and clustered.
                        note the lrms column is the Calpha RMSD of each fit from the first fit in its class
        
        
        """
        blurrer = StructureBlurrer()
        
        scorer = ScoringFunctions()
        
        cluster=Cluster()

        count=0
        dict_ensembl={}
        list_to_order=[]
        #print targetMap
        if targetMap==False:
            #targetMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print "WARNING:Need target map"
            sys.exit()
        if score not in ['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD']:
            print 'Incorrect Scoring Function: %s' % score
            print 'Please select from one of the following scoring functions: %s' % ', '.join(['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD'])
            sys.exit()
        
        
        targetMap=targetMap.copy()
        if score=='CCC':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]

                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None: score_mod=scorer.CCC_map(sim_map,targetMap,0.5*sim_map.fullMap.std(),cont_targetMap,2,True)[0]#CCC(sim_map,targetMap) 
                else: score_mod=scorer.CCC_map(sim_map,targetMap,0.0,0.0,True)[0]
                #else: score_mod=scorer.CCC(sim_map,targetMap)
                #'name_file','structure_instance','score','lrmsd','class'
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='LAP':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                score_mod=scorer.laplace_CCC(sim_map,targetMap) 
                #'name_file','structure_instance','score','lrmsd','class'
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='MI':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None: score_mod=scorer.MI(sim_map,targetMap,0.5*sim_map.fullMap.std(),cont_targetMap,1)
                else: score_mod=scorer.MI(sim_map,targetMap)
                list_to_order.append([name_mod,mod,score_mod,0,0])            
        if score=='NV':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                #These two values should be calculated for the experimental map, and only
                #need to be calculated once, at the beginning
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None: score_mod=scorer.normal_vector_score(targetMap,sim_map, cont_targetMap-(0.1*targetMap.std()), cont_targetMap+(0.1*targetMap.std()),Filter=None)
                else:
                    min_thr=targetMap.get_primary_boundary(mod.get_prot_mass_from_atoms(), targetMap.min(), targetMap.max())        
                    points=targetMap.get_point_map(min_thr,percentage=0.2)
                
                    max_thr=targetMap.get_second_boundary(min_thr, points, min_thr, targetMap.max(),err_percent=1)
                    score_mod=scorer.normal_vector_score(targetMap,sim_map, min_thr, max_thr,Filter=None)
                score_mod = 1 - (score_mod/3.14)
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='NV_Sobel':

            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None: score_mod=scorer.normal_vector_score(targetMap,sim_map, cont_targetMap-(0.1*targetMap.std()), cont_targetMap+(0.1*targetMap.std()),Filter='Sobel')
                else:
                    min_thr=targetMap.get_primary_boundary(mod.get_prot_mass_from_atoms(), targetMap.min(), targetMap.max())                
                    points=targetMap.get_point_map(min_thr,percentage=0.2)
                    max_thr=targetMap.get_second_boundary(min_thr, points, min_thr, targetMap.max(),err_percent=1)
                    score_mod=scorer.normal_vector_score(targetMap,sim_map, min_thr, max_thr,Filter='Sobel')
                score_mod = 1 - (score_mod/3.14)
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='NV_Laplace':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None: score_mod=scorer.normal_vector_score(targetMap,sim_map, cont_targetMap-(0.1*targetMap.std()), cont_targetMap+(0.1*targetMap.std()),Filter='Laplace')
                else:
                    min_thr=targetMap.get_primary_boundary(mod.get_prot_mass_from_atoms(), targetMap.min(), targetMap.max())                
                    points=targetMap.get_point_map(min_thr,percentage=0.2)
                    max_thr=targetMap.get_second_boundary(min_thr, points, min_thr, targetMap.max(),err_percent=1)
                    score_mod=scorer.normal_vector_score(targetMap,sim_map, min_thr, max_thr,Filter='Laplace')
                score_mod = 1 - (score_mod/3.14)
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='ENV':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                min_thr=targetMap.get_primary_boundary(mod.get_prot_mass_from_atoms(), targetMap.min(), targetMap.max())
                score_mod=scorer.envelope_score(targetMap,min_thr,mod)
                 #'name_file','structure_instance','score','lrmsd','class'
                list_to_order.append([name_mod,mod,score_mod,0,0])
        if score=='CD':        
            for mod1 in ensemble_list:
                count+=1
                name_mod=mod1[0]
                mod=mod1[1]
                sim_map = blurrer.gaussian_blur(mod, res_target_map,densMap=targetMap,sigma_coeff=sigma_coeff)
                if not cont_targetMap is None:
                    score_mod=scorer._surface_distance_score(sim_map,targetMap,0.5*sim_map.fullMap.std(),cont_targetMap,'Minimum')
                else:
                    min_thr=targetMap.get_primary_boundary(mod.get_prot_mass_from_atoms(), targetMap.min(), targetMap.max())
                    points=targetMap.get_point_map(min_thr,percentage=0.2)
                    max_thr=targetMap.get_second_boundary(min_thr, points, min_thr, targetMap.max(),err_percent=1)
                    score_mod=scorer.chamfer_distance(sim_map,targetMap, min_thr, max_thr, kdtree=None)
                    score_mod = 1/score_mod
                list_to_order.append([name_mod,mod,score_mod,0,0])

        if score in ['NV','NV_Sobel','NV_Laplace']:
            list_ordered=sorted(list_to_order, key=lambda x: x[2],reverse=True)#was false when NV was negative
        else:
            list_ordered=sorted(list_to_order, key=lambda x: x[2],reverse=True)
        if number_top_mod==0:
            if write==True:
                return cluster._print_results_cluster2(list_ordered,write)
            return list_ordered
        else:
            x=int(number_top_mod)
            if write==True:
                return cluster._print_results_cluster2(list_ordered[:x],write)
            return list_ordered[:x]
