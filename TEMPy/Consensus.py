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
from TEMPy.Cluster import  Cluster
from numpy import zeros,mean,median,asarray
from scipy.stats import mode
import sys
from collections import defaultdict

class Consensus:
    """A class to clustering an ensemble of structure instance"""
    
    def __init__(self):
        pass
        
    def _makedict_value(self,rankCCC):
        """
        private function used in Consensus Module.
        """
        #print rankCCC
        rank_dict={}
        for r in rankCCC:
            rank_dict[r[0]]=r[2]
        return rank_dict

    def _makedict(self,rank_score):
        """
        private function used in Consensus Module.
        """
		
        namerank_score=[mod[0] for mod in rank_score]
        d_rank={i:j for i,j in enumerate(namerank_score,start=1)}
        return d_rank

    def _makedict_list(self,list_score):
        """
        private function used in Consensus Module.
        """
        #print enumerate(rankCCC)
        d_rank={i:j for i,j in list_score}
        return d_rank

 
    def _printdict(self,dict_score):
        """
        private function used in Consensus Module.
        """
        for k,v in dict_score.items():
            print k,v

    def _modes(self,values):
        """
        private function used in Consensus Module.
        """
        count = defaultdict(int)
        for v in values:
            count[v] +=1
        best = max(count.values())
        print [k for k,v in count.items() if v == best]

    def _mode_here(self,arr):
        """
        private function used in Consensus Module.
        """
        m = max([arr.count(a) for a in arr])
        print [x for x in arr if arr.count(x) == m][0] if m>1 else None

           
    def vote_mode(self,ensemble_list,score_list,res_target_map,sigma_coeff,number_top_mod=0,write=False,targetMap=False):
        """
          
             Mode consensus scoring calculation between multiple "fits" using a user defined set of scores.
                           
                Arguments:
                    *ensemble_list*
                        Input list of Structure Instances.
                    *score_list*
                    	Input list of scoring function to use.
                        
                        See ScoringFunctions class for a list of the available Scoring Function.
                        E.g. set score='CCC' to use the Cross-correlation coefficient.
                        
                        Score option are:
                        
                        i    'CCC' - Cross-correlation coefficient; 
                        
                        ii    'LAP' - Laplacian-filtered cross-correlation coefficient:  useful for maps with resolutions worse than 10-15 A;
                        
                        iii   'MI' - Mutual information score: a good and robust score but relatively slow to calculate; 
                        
                        iv    'ENV' - Envelope score: the fastest score to calculate due to binarisation of the map. 
                        
                        v-vii  'NV','NV_Sobel','NV_Laplace'- Normal vector score: a vector-based surface superimposition score with or without Sobel/Laplace filter.

                        viii 'CD' - Chamfer Distance: a score used in computer vision algorithms as a fast similarity metric 
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
                    *targetMap*
                        Target Map Instance.

        """

        cluster=Cluster()
        list_dict=[]
        if targetMap==False:
            #targetMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print "WARNING:Need target map"
            sys.exit()
        score_select=[]
        for score in score_list:
            #check if score chosen are correct
            if score not in ['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD']:
                print 'Incorrect Scoring Function: %s' % score
                print 'Please select from one of the following scoring functions: %s' % ', '.join(['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD'])
                sys.exit()
            if score not in score_select:
                score_select.append(score)
            else:
                print 'Chose the %s twice' % score
                sys.exit()
        for score in score_list:
            print "******",score
            if score=='CCC':
                rankCCC=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())
                dictCCC=Consensus()._makedict(rankCCC)
                list_dict.append(dictCCC)
                Consensus()._printdict(dictCCC)
            elif score=='LAP':        
                rankLAP=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictLAP=Consensus()._makedict(rankLAP)
                list_dict.append(dictLAP)
                Consensus()._printdict(dictLAP)                        
            elif score=='MI':        
                rankMI=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictMI=Consensus()._makedict(rankMI)
                list_dict.append(dictMI)
                Consensus()._printdict(dictMI)
            
            elif score=='NV':        
                rankNV=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNV=Consensus()._makedict(rankNV)
                list_dict.append(dictNV)
                Consensus()._printdict(dictNV)
                
            elif score=='NV_Sobel':        
                rankNVS=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNVS=Consensus()._makedict(rankNVS)
                list_dict.append(dictNVS)
                Consensus()._printdict(dictNVS)               
            elif score=='NV_Laplace':        
                rankNVL=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNVL=Consensus()._makedict(rankNVL)
                list_dict.append(dictNVL)
                Consensus()._printdict(dictNVL)                  
            elif score=='ENV':        
                rankENV=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictENV=Consensus()._makedict(rankENV)
                list_dict.append(dictENV)           
                Consensus()._printdict(dictENV)  
            if score=='CD':        
                rankCD=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictCD=Consensus()._makedict(rankCD)
                list_dict.append(dictCD)           
                Consensus()._printdict(dictCD)         
        dict_count={}
        mxcinsensus = zeros(shape=(7,number_top_mod))
        for k,v in list_dict[0].items():
            dict_count[v]=[]
        for k in dict_count:
            for num in range(len(list_dict)):
                for k2,v2  in list_dict[num].items():
                    if k == v2:
                        dict_count[k].append(k2)
        dict_out={}
        for k,v in dict_count.items():
            median_list=median(v)
            m = max([v.count(a) for a in v])
            if m>1:
                mode_list=[x for x in v if v.count(x) == m][0]
                dict_out[k]=[median_list,mode_list]
            else:
                pass
            mode_list=max(set(v), key=v.count)
            
        sorted_dict = sorted(dict_out.items(), key=lambda x: x[1])
        print "**************"
        print "Consensus rank"
        for fit in sorted_dict:
            print fit[1],fit[0]
        return sorted_dict
    
    
    def _borda_score(self,list_rank,candidate,voters):
        """
        private function used in vote function.
    	It calculates the Borda count is a single-winner election method in which voters rank candidates in order of preference.
        """
        score=0
        for r in list_rank:
            score+=(candidate-r)*voters
        return score



    def vote(self,ensemble_list,score_list,res_target_map,sigma_coeff,number_top_mod=0,write=False,targetMap=False):
        """
          
             Borda consensus scoring calculation between multiple "fits" using a user defined set of scores.
             The Borda count is a single-winner election method in which voters rank candidates in order of preference.
                           
                Arguments:
                    *ensemble_list*
                        Input list of Structure Instances.
                    *score_list*
                    	Input list of scoring function to use.
                        
                        See ScoringFunctions class for a list of the available Scoring Function.
                        E.g. set score='CCC' to use the Cross-correlation coefficient.
                        
                        Score option are:
                        
                        i    'CCC' - Cross-correlation coefficient; 
                        
                        ii    'LAP' - Laplacian-filtered cross-correlation coefficient:  useful for maps with resolutions worse than 10-15 A;
                        
                        iii   'MI' - Mutual information score: a good and robust score but relatively slow to calculate; 
                        
                        iv    'ENV' - Envelope score: the fastest score to calculate due to binarisation of the map. 
                        
                        v-vii  'NV','NV_Sobel','NV_Laplace'- Normal vector score: a vector-based surface superimposition score with or without Sobel/Laplace filter.

                        viii 'CD' - Chamfer Distance: a score used in computer vision algorithms as a fast similarity metric 
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
                    *targetMap*
                        Target Map Instance.

        """
 

        cluster=Cluster()
        list_dict=[]
        candidate=len(ensemble_list)
        voters=len(score_list)
        if targetMap==False:
            #targetMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print "WARNING:Need target map"
            sys.exit()
        score_select=[]
        for score in score_list:
            #check if score chosen are correct
            if score not in ['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD']:
                print 'Incorrect Scoring Function: %s' % score
                print 'Please select from one of the following scoring functions: %s' % ', '.join(['CCC','LAP','MI','NV','NV_Sobel','NV_Laplace','ENV','CD'])
                sys.exit()
            if score not in score_select:
                score_select.append(score)
            else:
                print 'Chose the %s twice' % score
                sys.exit()
        for score in score_list:
            print "******",score
            if score=='CCC':
                rankCCC=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())
                #print rankCCC
                dictCCC=Consensus()._makedict(rankCCC)
                list_dict.append(dictCCC)
                Consensus()._printdict(dictCCC)
            elif score=='LAP':        
                rankLAP=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictLAP=Consensus()._makedict(rankLAP)
                list_dict.append(dictLAP)
                Consensus()._printdict(dictLAP)                        
            elif score=='MI':        
                rankMI=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictMI=Consensus()._makedict(rankMI)
                list_dict.append(dictMI)
                Consensus()._printdict(dictMI)
            elif score=='NV':        
                rankNV=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNV=Consensus()._makedict(rankNV)
                list_dict.append(dictNV)
                Consensus()._printdict(dictNV)
                for i in rankNV:
                	print i[0],i[2]
                
            elif score=='NV_Sobel':        
                rankNVS=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNVS=Consensus()._makedict(rankNVS)
                list_dict.append(dictNVS)
                Consensus()._printdict(dictNVS)               
            elif score=='NV_Laplace':        
                rankNVL=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictNVL=Consensus()._makedict(rankNVL)
                list_dict.append(dictNVL)
                Consensus()._printdict(dictNVL)                  
            elif score=='ENV':        
                rankENV=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictENV=Consensus()._makedict(rankENV)
                list_dict.append(dictENV)           
                Consensus()._printdict(dictENV)  
            if score=='CD':        
                rankCD=cluster.rank_fit_ensemble(ensemble_list,score,res_target_map,sigma_coeff,number_top_mod=number_top_mod,targetMap=targetMap.copy())   
                dictCD=Consensus()._makedict(rankCD)
                list_dict.append(dictCD)           
                Consensus()._printdict(dictCD)         
        dict_count={}
        #dict with [0,0,0,0,0,0,0,score_mod] possibility CCC,MI ... so keep order. add a colum for 'goal done -goal concived'
        #in our case how many time is 1st.
        #for d_rank in list_dict:
                #sorted_dict = sorted(dict.items(), key=lambda x: x[1])
                #print sorted_dict
        mxcinsensus = zeros(shape=(7,number_top_mod))
        for k,v in list_dict[0].items():
            dict_count[v]=[]
        for k in dict_count:
            #print 'k',k
            for num in range(len(list_dict)):
                for k2,v2  in list_dict[num].items():
                    if k == v2:
                        dict_count[k].append(k2)
        dict_out={}
        for k,v in dict_count.items():
           # print k
           # print sum(v)
            #print mean(v)
            median_list=median(v)
            #print mode(v)
            #most_frequent=mode(v)[0][0]
            #print 'mode most freq',most_frequent
            #Consensus().modes(v)
            #Consensus().mode_here(v)
            #print v
            borda_score=Consensus()._borda_score(v,candidate,voters)
            dict_out[k]=[borda_score,v]
            
        sorted_dict = sorted(dict_out.items(), key=lambda x: x[1][0],reverse=True)
        print "**************"
        print "Consensus rank"
        line=''
        line+="Borda_score\t"
        for score in score_list:
            line+='%s\t'%score
        line+="Fit\n"
        count=0
        for fit in sorted_dict:
            count+=1
            line+='%s\t'%count
            b=fit[1][0]
            line+='%s\t'%b
            for s in fit[1][1]:
                line+='%s\t'%s
            m=fit[0]
            line+='%s\n'%m
        print line
        return sorted_dict

#need to make it more elegant this come from private scripting.
    def vote_list(self,score_lists):
        """
          
             Borda consensus scoring calculation between multiple "fits" using a user defined set of scores.
             The Borda count is a single-winner election method in which voters rank candidates in order of preference.
                           
                Arguments:
                    *ensemble_list*
                        Input list of Structure Instances.
                    *score_list*
                    	Input list of list. Each list is a list of Structure Instances associated with a score.
        """
        
        dict_count={}
        list_dict=[]
        candidate=[]
        voters=len(score_lists)
        for i in score_lists:
        	candidate.append(len(i))
        for list_score in  score_lists:
            dictScore=Consensus()._makedict(list_score)
            list_dict.append(dictScore)
        for k,v in list_dict[0].items():
            dict_count[v]=[]
        for k in dict_count:
            #print 'k',k
            for num in range(len(list_dict)):
                for k2,v2  in list_dict[num].items():
                    if k == v2:
                        dict_count[k].append(k2)
        dict_out={}
        for k,v in dict_count.items():
            #print v
            #v = asarray(v)
            #print v
            #median_list=median(v)
            borda_score=Consensus()._borda_score(v,candidate[0],voters)
            dict_out[k]=[borda_score]
            
        sorted_dict = sorted(dict_out.items(), key=lambda x: x[1][0],reverse=True)
        print "**************"
        print "Consensus rank"
        line=''
        line+="Borda_score\t"
        count=0
        for score in score_lists:
        	count+=1
        	line+='%s\t'%count
        line+="Fit\n"
        count=0
        for fit in sorted_dict:
            count+=1
            line+='%s\t'%count
            b=fit[1][0]
            line+='%s\t'%b
            for s in fit[1][1]:
                line+='%s\t'%s
            m=fit[0]
            line+='%s\n'%m
        print line
        return sorted_dict


#     def vote_list(self,score_lists):
# 
#         dict_count={}
#         list_dict=[]
#         for list_score in  score_lists:
#             dictScore=Consensus()._makedict_list(list_score)
#             list_dict.append(dictScore)
#         for k,v in list_dict[0].items():
#             dict_count[v]=[]
#         for k in dict_count:
#             for num in range(len(list_dict)):
#                 for k2,v2  in list_dict[num].items():
#                     if k == v2:
#                         dict_count[k].append(k2)
#         dict_out={}
#         for k,v in dict_count.items():
#            # print k
#            # print sum(v)
#             #print mean(v)
#             #print median(v)
#             most_frequent=mode(v)[0][0]
#             dict_out[k]=most_frequent
#         sorted_dict = sorted(dict_out.items(), key=lambda x: x[1])
#         print "**************"
#         print "Consensus rank"
#         for fit in sorted_dict:
#             print fit[1],fit[0]
#         return sorted_dict
# 
            
            
            
            
            
            