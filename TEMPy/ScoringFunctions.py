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

from TEMPy.StructureBlurrer import StructureBlurrer
import math
from numpy import sum as numsum, copy as npcopy,mean as npmean, log10 as np_log10
from numpy import square,sqrt,absolute,histogram,argwhere,amin,count_nonzero,shape,size, array as nparray,\
                  transpose, mgrid,indices,meshgrid,nonzero,real,searchsorted,newaxis,where,matrix,ravel,ma,\
                  amax,ones,arange,floor,ceil,zeros, conjugate
from scipy.ndimage.interpolation import map_coordinates,spline_filter
from scipy.fftpack import fftn, ifftn, fftshift, fftfreq, ifftshift
#from scipy import weave
#from scipy.weave import converters
from scipy.spatial import KDTree
import sys
import itertools

class ScoringFunctions:
    """ 
    
    A class implementing various scoring functions used in density fitting. 
    Reference:
    Vasishtan and Topf (2011) Scoring functions for cryoEM density fitting.
    J Struct Biol 174:333-343.
    
    """
    def __init__(self):
        pass


    def _overlap_map_samebox(self,map1,map2):
        """
        
       volume overlap within 2 maps with same box size

        Return:
           % of overlap       
        
        """
        
        b=map1.fullMap
        binmap1=map1.fullMap>0.0
        binmap2=map2.fullMap>0.0
        mask_array=(binmap1*binmap2)>0.0
        return[count_nonzero(binmap1),count_nonzero(binmap2),count_nonzero(mask_array),mask_array.size]


    def _overlap_map_array(self,map_target,map_target_threshold,map_probe,map_probe_threshold):
        """
            mask maps with 2 cut-off map_target_threshold and map_probe_threshold (vol thr.)
            
            return:
            mask array where both are true.
            
        """
        
        binmap1=map_target.fullMap>float(map_target_threshold)
        binmap2=map_probe.fullMap>float(map_probe_threshold)
        mask_array=(binmap1*binmap2)>0
        return mask_array
    
    #add by AJP
    def calculate_map_threshold(self,map_target):
        try:
            peak,ave,sigma = map_target._peak_density()
            vol_threshold = float(ave)+(2.0*float(sigma))
        except:
            if len(map_target.header)==0:
                #amin = map_target.min()
                #amax = map_target.max()
                amean = map_target.mean()
                rms = map_target.std()
                vol_threshold = float(amean)+(1.5*float(rms))
            else:
                #amin = map.header[19]
                #amax = map.header[20]
                amean = map_target.mean()
                rms = map_target.std()
                vol_threshold = float(amean)+(1.5*float(rms))
        
        return vol_threshold
        

    def mapComparison(self, map_target, map_probe):
        """
        
        Compare the properties (sampling rate, box size and origin) of two maps 
        Arguments:
            *map_target, map_probe*
                Map instances to compare.
        Return:
            True if the map properties are the same between two maps, False otherwise.
        
        """
        
        if (map_target.apix - map_probe.apix < 1E-6) and map_target.box_size() == map_probe.box_size(): 
            if round(map_target.origin[0],2) == round(map_probe.origin[0],2) and round(map_target.origin[1],2) == round(map_probe.origin[1],2) and round(map_target.origin[2],2) == round(map_probe.origin[2],2):
                return True
            else:
                return False
        else: return False

    def _failed_match(self):
        print("Warning: can't match the map at the moment, use map with same box size.") #comment all out!
        sys.exit()
    
    def _CCC_calc(self,m1,m2):
        arr1 = m1.view(float)
        arr2 = m2.view(float)
        nd = len(arr1.shape)
        if nd == 2 and len(arr1.shape)[1] == 0:
            nd = 1
        l = 1
        dim = zeros(3,dtype=int)
        for i in range(nd):
            l *= arr1.shape[i]
            dim[i] = arr1.shape[i]
        #l = len(arr1)
        corr = 0.0
        #dims = nparray(ltmp,dtype=int)
        code = """
        int k,j,i;
        float numer=0.0, var1=0.0, var2 = 0.0;
        if (nd == 1){ 
         for (int z=0; z<dim[0]; z++) {
          numer += arr1[z]*arr2[z];
          var1 += pow(arr1[z],2);
          var2 += pow(arr2[z],2); }
        }
        else if (nd == 3){
          for (int z=0; z<dim[0]; z++) {
            for (int y=0; y<dim[1]; y++) {
              for (int x=0; x<dim[2]; x++) {
                numer += ARR13(z,y,x)*ARR23(z,y,x);
                var1 += pow(ARR13(z,y,x),2);
                var2 += pow(ARR23(z,y,x),2);
              }
            }
          }
        }
        corr = (float) numer/sqrt(var1*var2);
        return_val = corr;      
        """

        # check
        # BEN - commented out
        #try:
        #    #print datetime.now().time()
        #    corr = weave.inline(code,['arr1','arr2','corr','nd','dim'],headers=["<math.h>"],verbose=0)
        #    #print datetime.now().time()
        #    corr = min(1.0,corr)
        #    corr = max(-1.0,corr)
        #    return corr
        #except:
        #    #print 'C++ scoring run failed!'

        #    return None
        return None

    # Cross correlation coefficient for the overlap (3), contoured (2) or complete map (1), added by APJ
    def CCC_map(self, map_target,map_probe,map_target_threshold=0.0,map_probe_threshold=0.0,mode=1,meanDist=False,cmode=True):
        """
        Calculate cross-correlation between two Map instances, for the overlap (3), contoured (2) or complete map (1).

        Arguments:
                *map_target, map_probe*
                        EMMap instances to compare.
                *map_target_threshold,map_probe_threshold*
                        EMMap threshold
                        if not given, use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold
                *mode*
                        3. calculation on the mask
                        2. calculation on contoured maps
                        1. calculation on complete map
                *meanDist*
                        True if the deviation from mean needs to be calculated
        """
        if self.mapComparison(map_target, map_probe):
            
            if not mode == 1:
                # calculate threshold if not given : 2* sigma can be used for experimental maps and 1*sigma for simulated?
                if map_target_threshold==0 and map_probe_threshold==0:
                    map_target_threshold=self.calculate_map_threshold(map_target)
                    map_probe_threshold=self.calculate_map_threshold(map_probe)
                # calculate contour overlap
                # contour the first map
                bin_map1 = map_target.fullMap > float(map_target_threshold)
                bin_map2 = map_probe.fullMap > float(map_probe_threshold)
                # percent calculated on the smaller contoured volume (can be changed)
                minim = numsum(bin_map1)
                minim2 = numsum(bin_map2)
                if minim2 < minim: minim = minim2
                mask_array = (bin_map1*bin_map2) > 0
                #print '>>', numsum(bin_map1),numsum(bin_map2),numsum(mask_array),minim
                if not minim == 0.0:perc_ovr = float(numsum(mask_array))/minim
                else:
                    perc_ovr = 0.0
                    print('No map overlap (Cross correlation score), exiting score calculation..')
                    return -1.0, 0.0
                if perc_ovr < 0.02: return -1.0, 0.0
            else: perc_ovr = 1.0

            # calculate CCC within volume of overlap
            if mode == 3:
                #mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
                if numsum(mask_array) == 0:
                    print ('No map overlap (Cross correlation score), exiting score calculation..')
                    return -1.0, 0.0
                map1_mask = map_target.fullMap[mask_array]
                map2_mask = map_probe.fullMap[mask_array]
                if meanDist:
                    map1_mask = map1_mask - npmean(map1_mask)
                    map2_mask = map2_mask - npmean(map2_mask)
                if cmode:
                    corr = self._CCC_calc(map1_mask.flatten(),map2_mask.flatten())
                    #print corr, numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask))), numsum(map1_mask * map2_mask)
                else: corr = None
                if corr is None:
                    return numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask))), perc_ovr
                else: return corr, perc_ovr
            # calculate CCC for contoured maps based on threshold
            elif mode == 2:
                #bin_map1 = map_target.fullMap > float(map_target_threshold)
                #bin_map2 = map_probe.fullMap > float(map_probe_threshold)
                map1_mask = map_target.fullMap*bin_map1
                map2_mask = map_probe.fullMap*bin_map2
                if meanDist:
                    map1_mask = map1_mask - npmean(map_target.fullMap[bin_map1])
                    map2_mask = map2_mask - npmean(map_probe.fullMap[bin_map2])
                    map1_mask = map1_mask*bin_map1
                    map2_mask = map2_mask*bin_map2
                else:
                    map1_mask = map_target.fullMap*bin_map1
                    map2_mask = map_probe.fullMap*bin_map2
                if cmode: corr = self._CCC_calc(map1_mask,map2_mask)
                else: corr = None
                #print corr, numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask)))
                if corr is None:
                    return numsum(map1_mask * map2_mask)/sqrt(numsum(square(map1_mask))*numsum(square(map2_mask))), perc_ovr
                else:
                    return corr, perc_ovr
            # calculate on the complete map
            if meanDist:
                if cmode: corr = self._CCC_calc(map_target.fullMap-npmean(map_target.fullMap),map_probe.fullMap-npmean(map_probe.fullMap))
                else: corr = None
                #print corr,numsum((map_target.fullMap-npmean(map_target.fullMap)) * (map_probe.fullMap-npmean(map_probe.fullMap)))/(sqrt(numsum(square(map_target.fullMap-npmean(map_target.fullMap)))*numsum(square(map_probe.fullMap-npmean(map_probe.fullMap)))))
                if corr is None:
                    return numsum((map_target.fullMap-npmean(map_target.fullMap)) * (map_probe.fullMap-npmean(map_probe.fullMap)))/(sqrt(numsum(square(map_target.fullMap-npmean(map_target.fullMap)))*numsum(square(map_probe.fullMap-npmean(map_probe.fullMap))))), perc_ovr
                else: return corr, perc_ovr
            if cmode: corr = self._CCC_calc(map_target.fullMap,map_probe.fullMap)
            else: corr = None
            #print corr, numsum(map_target.fullMap * map_probe.fullMap)/sqrt(numsum(square(map_target.fullMap))*numsum(square(map_probe.fullMap))), numsum(map_target.fullMap * map_probe.fullMap)
            if corr is None:
                return numsum(map_target.fullMap * map_probe.fullMap)/sqrt(numsum(square(map_target.fullMap))*numsum(square(map_probe.fullMap))), perc_ovr
            else: return corr, perc_ovr
        else:
            print("@@@ Maps could not be matched")
            return -1., 0.



    def CCC(self, map_target, map_probe):
        """
        
        Calculate cross-correlation between two Map instances.
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
        Return:
            CCC score
        
        """
        
        if self.mapComparison(map_target, map_probe):
            return (map_target.normalise().getMap()*map_probe.normalise().getMap()).mean()
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)
            #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()
 
   #TODO: check and delete the following
    '''
        ### Correlation coefficient about mean for the overlap mask
    def CCC_local(self, map_target,map_probe,map_target_threshold=0,map_probe_threshold=0):
        """
        
        Calculate cross-correlation about mean between two Map instances, for the overlap region.

        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            mean CCC score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)
            mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
            map_target_mask = map_target.fullMap[mask_array]
            map_target_mask = map_target_mask - float(map_target_mask.sum()/len(map_target_mask))
            map_probe_mask = map_probe.fullMap[mask_array]
            map_probe_mask = map_probe_mask - float(map_probe_mask.sum()/len(map_probe_mask))
            return absolute((map_target_mask * map_probe_mask)).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
            #return (map_target_mask * map_probe_mask).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
        else:
            self._failed_match()
                        #m1,m2 = self.matchMaps(map_target, map_probe)
                        #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()

        # MAIN: Cross correlation coefficient for the overlap (3), contoured (2) or complete map (1)
    
    def CCC_mask_zero(self, map_target,map_probe,map_target_threshold=0,map_probe_threshold=0):
        """
        
        Calculate cross-correlation about zero for the overlap region between two Map instances.
                                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            mean CCC score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)
            mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
            map_target_mask = map_target.fullMap[mask_array]
            map_probe_mask = map_probe.fullMap[mask_array]
            return (map_target_mask * map_probe_mask).sum()/sqrt(square(map_target_mask).sum()*square(map_probe_mask).sum())
        else:
            self._failed_match()
                        #m1,m2 = self.matchMaps(map_target, map_probe)
                        #return (m1.normalise().getMap()*m2.normalise().getMap()).mean()

    '''
    
    def LSF(self, map_target, map_probe):
        """
        
        Calculate least-squares between two Map instances.
        
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
        Return:
            least-squares value
            
        """
        if self.mapComparison(map_target, map_probe):
        	map_target, map_probe = map_target, map_probe
            
        else:
            self._failed_match()
        return ((map_target.getMap()-map_probe.getMap())**2).mean()

    def laplace_CCC(self, map_target, map_probe, prefil=(False, False)):
        """
        
        Calculate Laplacian cross-correlation between two Map instances.
        Based on (Chacon and Wriggers, 2002).
                
        Arguments:
            *map_target, map_probe*
                Map instances to compare.
            *prefil*
                2-tuple of boolean values, one for each map respectively.
                True if Map instance is already Laplacian-filtered. False otherwise.
        Return:
            Laplacian cross-correlation score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)

        if not prefil[0]:
            map_target = map_target.laplace_filtered()
        if not prefil[1]:
            map_probe = map_probe.laplace_filtered()
        map_target = map_target.normalise()
        map_probe = map_probe.normalise()
        return self.CCC(map_target, map_probe)

        # MAIN: normal vector score calculated on surface voxels derived by different methods
    def normal_vector_score(self, map_target, map_probe, primary_boundary, secondary_boundary=0.0,Filter=None):
        """
        
        Calculate the Normal Vector Score between two Map surfaces.
        Based on 3SOM algorithm (Ceulemans and Russell, 2004). 
        
                
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare. map_target is the target map.
            *primary_boundary, secondary_boundary*
                If a filter is selected, just input a contour level as primary threshold.
                Otherwise, need to run get_primary_boundary and get_second_boundary based on map target.
            *Filter*
            	Filter to use:  
            		
            		i Sobel Filter (Filter=='Sobel')
            		
            		ii Laplace Filter (Filter=='Laplace')
            		
            		iii Minimum Filter (Filter=='Minimum')
            		
            		iv Mean Filter (Filter=='Mean')
            		
            		
        Return:
            Normal vector score.
            
        """
        
        if Filter not in ['Sobel','Laplace','Mean','Minimum',None]:
        	print("Incorrect name of filter: " + Filter)
        	print("Select one of the following Filters if applicable: " + ''.join(['Sobel','Laplace']))
        	sys.exit()
        
        scores = []
        if not self.mapComparison(map_target, map_probe):
            #map_target, map_probe = self.matchMaps(map_target, map_probe)
            self._failed_match()
        assert isinstance(primary_boundary,float)
        assert isinstance(secondary_boundary,float)
        #print "fff", primary_boundary, secondary_boundary
        if primary_boundary > secondary_boundary:
            temp_thr = secondary_boundary
            secondary_boundary = primary_boundary
            primary_boundary = temp_thr
                    
        points = argwhere((map_target.fullMap > primary_boundary) & (map_target.fullMap < secondary_boundary))
        if Filter=='Sobel':
        	# sobel filter surface
            map1_surface = map_target._sobel_filter_contour(primary_boundary)
            points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
        elif Filter=='Laplace':
            # sobel filter surface
            map1_surface = map_target._laplace_filtered_contour(primary_boundary)
            points = argwhere(map1_surface.fullMap > (map1_surface.max()/2.0))
        elif  Filter=='Minimum':
            # the filter returns points touching surface (zeros)
            #map1_surface = map_target._surface_minimum_filter(float(primary_boundary))
            map1_surface = map_target._surface_minimum_filter(float(primary_boundary))
            points = argwhere(map1_surface == 1)
        elif Filter=='Mean':
            # the filter returns points from protrusions/curved surfaces
            map1_filter = map_target._surface_features(float(primary_boundary))
            # to extract points with filtered values less than a cut-off
            # more finer the bins are, more precise will be number of points chosen; not very crucial
            bin_test = [0.0001]
            for ii in range(1,41): bin_test.append(0.025*ii)
            freq_test = histogram(map1_filter.fullMap,bin_test)[0]
            sum_freq = 0.0 
            for fr in range(len(freq_test)):
                sum_freq += float(freq_test[fr])
                if sum_freq/numsum(freq_test) > 0.05 and bin_test[fr+1] >= 0.3:
                    t1 = bin_test[fr+1]
                    break
                if sum_freq/numsum(freq_test) > 0.10 or sum_freq > 100000:
                    t1 = bin_test[fr+1]
                    break
            points = argwhere((map1_filter.fullMap > 0.0) & (map1_filter.fullMap < t1))
        #C++ calculation
        flagc = 1
        try:    
            vecnorm_target =  map_target._get_normal_vector(points)
            vecnorm_probe = map_probe._get_normal_vector(points)
        except:
            flagc = 0
        if vecnorm_target is None or vecnorm_probe is None: flagc = 0
        ct = 0
        
        if flagc == 1:
            for l in range(len(vecnorm_target)):
                ct += 1
                nvec = vecnorm_target[l]
                ovec = vecnorm_probe[l]
                ### add max value for regions of null variation
                if (nvec[0] == 0. and nvec[1] == 0. and nvec[2] == 0.):
                    if (ovec[0] == 0. and ovec[1] == 0. and ovec[2] == 0.0):
                        continue
                    else:
                        scores.append(3.14)
                        continue
                else:
                    if (ovec[0] == 0. and ovec[1] == 0. and ovec[2] == 0.):
                        scores.append(3.14)
                        continue
                
                try:
                    dotprod = ovec[0] * nvec[0] + ovec[1] * nvec[1] + ovec[2] * nvec[2]
                    den = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2) * sqrt(ovec[0]**2 + ovec[1]**2 + ovec[2]**2)
                    if abs(dotprod-den) < 0.00001:
                        ang = 0.0
                    else:
                        ang = math.acos(min(max(dotprod/den,-1.0),1.0))
                    if den == 0.0: print(dotprod, den, nvec, ovec)
                    scores.append(abs(ang))   
                except ValueError:
                    print('Error: Angle could not be calculated: ', nvec,' ', ovec)
            #print scores[-10:]
            if len(scores) == 0:
                print("There are no points to be scored! The threshold values or the number of points to be considered needs to be changed.")
                return None
            else:
                if sum(scores) == 0:
                    return 0.0
                else:
                    #return 1-(sum(scores)/(len(points)*3.14)) #in this way go from 1 to 0
                    return 1-(sum(scores)/(len(points)*3.14))
        scores = []
        ct1 = 0
        if flagc == 0:    
            for v in points:
                n_vec = map_target.get_normal_vector(v[2],v[1],v[0])
                o_vec = map_probe.get_normal_vector(v[2],v[1],v[0])
                ct1 += 1
                ### add max value for regions of null variation
                if (n_vec.x == -9 and n_vec.y == -9 and n_vec.z == -9):
                    if (o_vec.x == -9 and o_vec.y == -9 and o_vec.z == -9):
                        continue
                    else:
                        scores.append(3.14)
                        continue
                else:
                    if (o_vec.x == -9 and o_vec.y == -9 and o_vec.z == -9):
                        scores.append(3.14)
                        continue
                try:
                    scores.append(abs(n_vec.arg(o_vec)))
                except ValueError:
                    print('Error: Angle between '+ str(n_vec) +', '+ str(o_vec) +' for point %d, %d, %d cannot be calculated.',v.x,v.y,v.z)
        if len(scores) == 0:
            print ("There are no points to be scored! The threshold values or the number of points to be considered needs to be changed.")
        else:
            if sum(scores) == 0:
                return 0
            else:
                #return 1-(sum(scores)/(len(points)*3.14)) #in this way go from 1 to 0
                return 1-(sum(scores)/(len(points)*3.14))


    def get_partial_DLSF(self, num_of_points, map_target, map_probe):
        """
        
        Calculate the DLSF score between two Map instances.
        The DLSF is similar to the LSF; 
        whereas the LSF compares absolute density values, 
        the DLSF compares the difference between pairs of values. 
    
        Arguments:
            *map_target, map_probe*
                the two Map instances to compare.
            *num_of_points*
                number of significant points.
        Return:
            DLSF score        
        """
        
        if not self.mapComparison(map_target, map_probe):
            #map_target, map_probe = self.matchMaps(map_target, map_probe)
            return "can't Match the map"
        #print "fff", primary_boundary, secondary_boundary

        map_target_sig_pairs=map_target._get_random_significant_pairs(int(num_of_points))
        otherMap=map_probe
        score = 0.0
        for p in map_target_sig_pairs:
            z1 = p[0]
            y1 = p[1]
            x1 = p[2]
            z2 = p[3]
            y2 = p[4]
            x2 = p[5]
            dens = p[6]
            prot_dens = otherMap.fullMap[z1][y1][x1] - otherMap.fullMap[z2][y2][x2]
            score += (dens-prot_dens)**2
        return score/map_target.fullMap.size
        
    def _MI(self, map_target, map_probe, layers=20):
        """ 
        
        Calculate the mutual information score between two Map instances.
                     
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *layers*
                Number of layers used to bin the map. Default is 20  as in Shatsky et al., 2008.
           Return:
            MI score
        
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = self.matchMaps(map_target, map_probe)
        score = 0
        m1_levels = (m1.max()-m1.min())/layers
        m2_levels = (m2.max()-m2.min())/layers
        for x in range(layers):
            for y in range(layers):
                m1_level_map = (m1.getMap() >= m1.min()+(x*m1_levels))*(m1.getMap() <= m1.min()+((x+1)*m1_levels))
                m2_level_map = (m2.getMap() >= m2.min()+(y*m2_levels))*(m2.getMap() <= m2.min()+((y+1)*m2_levels))
                comb_level_map = m1_level_map*m2_level_map
                p_m1 = float(m1_level_map.sum())/m1_level_map.size
                p_m2 = float(m2_level_map.sum())/m2_level_map.size
                p_comb = float(comb_level_map.sum())/comb_level_map.size
                if p_comb == 0:
                    mi_score = 0.0
                else:
                    #print p_comb, p_m1, p_m2, p_comb/(p_m1*p_m2), math.log(p_comb/(p_m1*p_m2),2)
                    mi_score = p_comb*math.log(p_comb/(p_m1*p_m2), 2)
                score += mi_score
        return score
    
    def _MI_C(self,m1,m2,layers1=20,layers2=20,N=0,lc1=0.0,lc2=0.0):
        #from datetime import datetime
        #print datetime.now().time()
        ly1 = int (layers1)
        ly2 = int (layers2)
        # input 3D arrays
        arr1 = (m1).view(float)
        arr2 = (m2).view(float)
        nz = int(arr1.shape[0])
        ny = int(arr1.shape[1])
        nx = int(arr1.shape[2])
        # min and max to set left and right bound
        ma1 = ma.masked_less_equal(arr1,lc1,copy=False)
        min1 = float(ma1.min())
        max1 = float(ma1.max())
        #print min1,max1,amin(m1[msk]),amax(m1[msk])
        #min1 = float(amin(m1[msk]))
        #max1 = amax(m1[msk])
        ma2 = ma.masked_less_equal(arr2,lc2,copy=False)
        min2 = float(ma2.min())
        max2 = float(ma2.max())
        #print min2,max2
        #min2 = float(amin(m2[msk]))
        #max2 = amax(m2[msk])
        min1 = float(min1-((max1-min1)/layers1)*0.0001)
        min2 = float(min2-((max2-min2)/layers2)*0.0001)
        # bin width
        step1 = (max1-min1)/float(layers1)
        step2 = (max2-min2)/float(layers2)
        # histogram freq in bins
        freq1 = zeros(layers1,dtype=float)
        freq2 = zeros(layers2,dtype=float)
        comb_freq = zeros((layers1,layers2),dtype=float)
        code = """
        int i,j,k,s1=0,s2=0;
        float p1=0.0, p2=0.0, pcomb = 0.0,Hxy=0.0,Hy=0.0,Hx=0.0;
        float va1,va2;
        /*long index = 0;
        long indexend = nz * ny * nx;
        while (index < indexend){
              va1  = arr1[index];
              va2 = arr2[index];*/
        /* use 3d array loop */
        for (int z=0; z<nz; z++) {
          for (int y=0; y<ny; y++) {
            for (int x=0; x<nx; x++) {
              va1 = ARR13(z,y,x);
              va2 = ARR23(z,y,x);
              for (i=0; i<ly1; i++) 
              {
                if ((va1 > (min1+ i*step1)) && (va1 <= (min1+(i+1)*step1))) 
                {
                  FREQ11(i) += 1.0;
                  s1 += 1;
                  break;
                 }
              }
              if (i == ly1) i = i-1;
              for (j=0; j<ly2; j++)
              {
                if ((va2 > (min2+j*step2)) && (va2 <= (min2+(j+1)*step2))) 
                {
                  FREQ21(j) += 1.0;
                  s2 += 1;
                  COMB_FREQ2(i,j) += 1.0;
                  break;
                }
              }
            /*index ++;*/
            }

          }
        }

        for (i=0; i<ly1; i++){
         p1 = FREQ11(i)/(float) s1;
         /*std::cout << s1 << ' ' << s2 << std::endl;*/
         for (j=0; j<ly2; j++){
           p2 = FREQ21(j)/(float) s2;
           pcomb = COMB_FREQ2(i,j)/(float) s1;
           if (pcomb != 0.0) Hxy += (-pcomb*log2(pcomb));
           if ((i == 0) && (p2 != 0.0)) Hy += (-p2*log2(p2));
         }
         if (p1 != 0.0) Hx += (-p1*log2(p1));
        }
        /*std::cout << Hxy << ' ' << Hx << ' ' << Hy << ' ' << std::endl;*/
        if (N == 1) {
          if (Hxy != 0.0) return_val = (Hx+Hy)/Hxy;
          else return_val = 0.0;
          }
        else return_val = Hx+Hy-Hxy;
            
        """


        # check
        # BEN Commented out due to weave
        #try:
            #print datetime.now().time()
        #    mi = weave.inline(code,['arr1','arr2','ly1','ly2','N','freq1','freq2','comb_freq','nz','ny','nx','step1','step2','min1','min2'],headers=["<math.h>"],verbose=0)
        #    #print datetime.now().time()
        #    mi = max(0.0,mi)
        #    return mi
        #except:
        #    #print 'C++ MI scoring run failed!'
        #    return None
        return None
    
        #Faster version of MI, in the overlap region (3) or complete density (1), added by APJ
    def MI(self, map_target, map_probe, map_target_threshold=0.0, map_probe_threshold=0.0, mode=1, layers1=None,layers2=None, weight=False,cmode=True):
        """ 
        
        Calculate the mutual information score between two Map instances.
                     
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold, map_probe_threshold*
                Thresholds used for contouring
            *mode*
                1. use complete map for calculation
                3. use overlap region for calculation
            *layers1, layers2*
                Number of layers used to bin the maps. Default is 20  as in Shatsky et al., 2008.
           Return:
            MI score
        
        """
        if not self.mapComparison(map_target, map_probe):
            #m1, m2 = map_target, map_probe
        #else:
            self._failed_match()
        # calculate threshold if not given : 2* sigma can be used for experimental maps and 1*sigma for simulated?
        if map_target_threshold==0.0:
            map_target_threshold=self.calculate_map_threshold(map_target)
        if map_probe_threshold==0.0:
            map_probe_threshold=self.calculate_map_threshold(map_probe)
        # calculation on the complete map
        if mode == 1:
            if weight: wt = 1
            else: wt = 0
            if layers1 is None:
                layers1 = 20
            if layers2 is None:
                layers2 = 20
            min1 = amin(map_target.fullMap) - 0.00001*(amax(map_target.fullMap)-amin(map_target.fullMap))
            min2 = amin(map_probe.fullMap) - 0.00001*(amax(map_probe.fullMap)-amin(map_probe.fullMap))
            if cmode: mic = self._MI_C(map_target.fullMap,map_probe.fullMap,layers1,layers2,wt,min1,min2)
            else: mic = None
            if not mic == None: return mic
            # digitize whole map based on layers
            map1_bin = map_target._map_digitize(map_target.min(),layers1,True)
            map2_bin = map_probe._map_digitize(map_probe.min(),layers2,True)
            bins1 = []
            for i in range(layers1+2): bins1.append(i)
            bins2 = []
            for i in range(layers2+2): bins2.append(i)
            # calculate frequency of bins
            map1_freq = histogram(map1_bin.fullMap,bins1)[0][1:]
            map2_freq = histogram(map2_bin.fullMap,bins2)[0][1:]
        elif mode == 3:
            # For score within masked region, the background is a bit ambiguous because low densities are overrepresented
            mask_array = self._overlap_map_array(map_target,map_target_threshold,map_probe,map_probe_threshold)
            if numsum(mask_array) == 0:
                print('No map overlap (Mutual information score), exiting score calculation..')
                return 0.0
            # sturges rule provides a way of calculating number of bins : 1+math.log(number of points)
            if layers1 is None:
                try: layers1=int(1+math.log(numsum(mask_array),2))
                except ValueError:
                    print('No map overlap (Mutual information score), exiting score calculation..')
                    return 0.0
            if layers2 is None:
                try: layers2=int(1+math.log(numsum(mask_array),2))
                except ValueError:
                    print('No map overlap (Mutual information score), exiting score calculation..')
                    return 0.0
            layers1 = max(layers1,15)
            layers2 = max(layers2,15)
            if weight: wt = 1
            else: wt = 0
            if cmode: mic = self._MI_C(nparray(map_target.fullMap*mask_array),nparray(map_probe.fullMap*mask_array),layers1,layers2,wt)
            else: mic = None
            if not mic == None: return mic
            # digitize masked map based on layers
            map1_bin = map_target.copy()
            map2_bin = map_probe.copy()
            map1_bin.fullMap = map1_bin.fullMap*mask_array
            map2_bin.fullMap = map2_bin.fullMap*mask_array
            map1_bin = map1_bin._map_digitize(map_target.fullMap[mask_array].min(),layers1,True)
            map2_bin = map2_bin._map_digitize(map_probe.fullMap[mask_array].min(),layers2,True)
            # make sure the outside region is filled with zeros
            map1_bin.fullMap = map1_bin.fullMap*mask_array
            map2_bin.fullMap = map2_bin.fullMap*mask_array
            #background frequencies from the whole map
            bins1 = []
            for i in range(layers1+2): bins1.append(i)
            bins2 = []
            for i in range(layers2+2): bins2.append(i)
            # calculate frequency of bins
            map1_freq = histogram(map1_bin.fullMap,bins1)[0][1:]
            map2_freq = histogram(map2_bin.fullMap,bins2)[0][1:]

        score = 0.0
        total = 0

        if numsum(map1_freq) == 0:
            print('No map overlap (Mutual information score), exiting score calculation..')
            return 0.0
        if numsum(map2_freq) == 0:
            print('No map overlap (Mutual information score), exiting score calculation..')
            return 0.0
        list_overlaps = []
        for x in range(layers1):
            mask_array = map1_bin.fullMap == float(x+1)
            overlap_freq =  histogram(map2_bin.fullMap[mask_array],bins2)[0][1:]
            total += float(numsum(overlap_freq))
            list_overlaps.append(overlap_freq)

        if total == 0:
            print('No map overlap (Mutual information score), exiting score calculation..')
            return 0.0
        enter = 0
        Hxy = 0.0
        Hx = 0.0
        Hy = 0.0
        mi_score = 0.0
        p_comb = 0.0
        #print numsum(map1_freq), numsum(map2_freq), total
        for x in range(layers1):
            # probability of occurrence of x
            p_m1 = map1_freq[x]/float(numsum(map1_freq))
            for y in range(layers2):
                enter = 1
                # probability for overlap of bins x and y
                p_comb = list_overlaps[x][y]/total
                # probability of occurrence of y
                p_m2 = map2_freq[y]/float(numsum(map2_freq))
                #if p_m1 == 0.0 or p_m2 == 0.0:
                #       mi_score = 0.0
                #       continue
                if p_comb == 0:
                    mi_score = 0.0

                else:
                    # p_m1 and p_m2 (background probabilties can be non-zero when p_comb=0), so the entropy based definition may be used
                    ## mi_score = p_comb*math.log(p_comb/(p_m1*p_m2), 2)
                    Hxy += -p_comb*math.log(p_comb, 2) # joined entropy
                score += mi_score
                if x == 0 and not p_m2 == 0.0: Hy += (-p_m2*math.log(p_m2, 2))
            if not p_m1 == 0.0: Hx += (-p_m1*math.log(p_m1, 2))
        if enter == 1:
            # normalised MI (Studholme et al.) is used to account for overlap of 'contours'
            # MI = Hx+Hy-Hxy & NMI = Hx+Hy/Hxy
            if weight:
                if Hxy == 0.0: return 0.0
                return (Hx+Hy)/Hxy
            return Hx+Hy-Hxy#score
        else: return None

        # MAIN: Faster version of MI, in the overlap region (3) or map contour (2) or complete density (1)
        
    def _hausdorff_list(self, primary_boundary, secondary_boundary, kdtree, map_probe):
        """
        
        This is for the chamdef distance def chamfer_distance, min max density value that define the surface of the protein
        
        Arguments:
        
            *kdtree* (there are 2 of them in numpy one Cbased on py-based, the latter is better, ctrl) this have to be one of the input.
                    kdtree from map_target 
            *primary_boundary, secondary_boundary*  need to run get_primary_boundary and get_second_boundary for map_probe
            
            NOTE: if you keep the kdtree as parametre out os less time consuming as building it takes time.
            
        """
        points = map_probe.get_pos(primary_boundary, secondary_boundary)
        #print "HERE POINTS",points
        return kdtree.query(points)[0] #kdtree give 2 list 0=distance 1=actual points
        

    def chamfer_distance(self, map_target, map_probe, primary_boundary, secondary_boundary, kdtree=None):
        """ 
        
        Calculate the chamfer distance Score between two Map instances. 
        NOT RACCOMANDED.      
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *primary_boundary*
                is the value returned by get_primary_boundary for map_probe
            *secondary_boundary*  
                is the value returned by get_second_boundary for map_probe
            *kdtree* 
                If set True it is possible to choose between the option of kdtree in numpy 
                The one that is py-based is a better choice.
        
        """
        
        if self.mapComparison(map_target, map_probe):
            m1, m2 = map_target, map_probe
        else:
            self._failed_match()
            #m1,m2 = matchMaps(map_target, map_probe)
        print("here")
        if kdtree:
            return self._hausdorff_list(primary_boundary, secondary_boundary, kdtree, m2).mean()
        else:
        	print(m1,primary_boundary, secondary_boundary)
        	kdtree = m1.makeKDTree(primary_boundary, secondary_boundary) #if you don't assine it wil be build one kdtree
        	if kdtree==None:
        		print("Error. No points selected, change boundary parameters.")
        		sys.exit()
      
        return self._hausdorff_list(primary_boundary, secondary_boundary, kdtree, m2).mean()#mean distance to the nearest neighbour 

        # CHAMFER DISTANCE SCORE based on a defined surface based on modes
    def _surface_distance_score(self,map_target,map_probe,map_target_threshold1=0.0,map_probe_threshold=0.0,Filter=None,map_target_threshold2=0.0,weight=False):
        """ 
        
        Calculate the chamfer distance Score between two Map instances. 
              
                        
        Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold1*
                contour threshold of the target map. 
                This value is used the primary boundary if map_target_threshold2 is given.
            *map_probe_threshold*
                contour threshold for the probe map.
            *Filter*
                definition of the surface:
                    1) None : surface defined by known boundaries - map_target_threshold1 & map_target_threshold2
                    If the boundaries are not known and target&probe map contour levels are known:
                        2) Std : to define the boundaries, contour level +- 5%sigma is calculated. 
                           5%sigma is used to limit the number of points picked as surface. 
                           For small maps, higher values (eg: 10%sigma) can be used.
                        3) Mean: a mean filter is applied on the binary contour mask over a long window. 
                           The resulting mask has values between 0 and 1. 
                           Points with values less than 0.3 is used to represent surface. 
                           As the average is calculated on a long window, highly exposed surface points \ 
                             have very low values and partially exposed surfaces/grooves have relatively higher values. 
                           This definition is useful especially when the map surface has many features/projections.
                        4) Minimum: a minimum filter is applied on a binary contour mask to locate surface points.
                           Voxels surrounded by points outside the contour (zeroes) are detected as surface.
                           Voxels surrounded by points outside the contour (zeroes) are detected as surface.
                        5) Sobel: sobel filter is applied on the map to detect high density gradients. 
                           Before applying the sobel filter, it is important to reduce the noise density \
                             and large variations (gradients) in the noise region.
            *weight*
                If set true, the distances between the surface points is normalized in a way similar to GDT (Zemla 2007)\
                  calculation for atomic co-ordinate alignments.

        """
        # check if both maps are on the same grid       
        if not self.mapComparison(map_target, map_probe):
            print("@@@ Maps could not be matched")
            return -999.
        # if the boundaries are known, calculate the kdtree
        if Filter == None:
            kdtree = map_target.makeKDTree(map_target_threshold1,map_target_threshold2)
            probe_points = map_probe.get_pos(map_target_threshold1, map_target_threshold2)

        # surface based on contour density thresholds for target and probe. 5% sigma is used to define boundaries.
        elif Filter == 'Std':
            # argwhere returns points as z,y,x, in the same way the map array dimensions are defined.
            target_points = argwhere((map_target.fullMap > (float(map_target_threshold1)-(map_target.std()*0.10))) & (map_target.fullMap < (float(map_target_threshold1)+(map_target.std()*0.10))))
            probe_points = argwhere((map_probe.fullMap > (float(map_probe_threshold)-(map_probe.std()*0.10))) & (map_probe.fullMap < (float(map_probe_threshold)+(map_probe.std()*0.10))))
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            if len(target_points) == 0 or len(probe_points) == 0:
              print ('Surface detection failed (Std filter), exiting..')
              return None
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return None
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return None
        elif Filter == 'Mean':
            map1_filter = map_target._surface_features(float(map_target_threshold1))
            map2_filter = map_probe._surface_features(float(map_probe_threshold))
            # define surface based on the filtered mask values.
            # points with values less than 0.3 are usually preferred. But in some cases like viruses, most surface points are highly exposed and \ 
            # a large number of points are returned and the calculation becomes slow.
            # Hence an additional filter is added: the maximum allowed points is 10% of box size. 
            # The minimum number of points is kept as 7%. This mode is less sensitive to the number of surface points chosen \
            # as the extent of exposure is used for defining surface. Hence thick surface is not usually required.

            # calculate frequencies in bins for filtered mask. 
            # The smaller the bins, more precise will be the calculation of points allowed based on percent of points chosen.
            # As this is just an additional filter and doesn't affect the calculations drastically, 40 bins are used to calculate frequencies.
            bin_test = [0.0001]
            for ii in range(1,41): bin_test.append(0.025*ii)
            freq_test = histogram(map1_filter.fullMap,bin_test)[0]
            map1_filled = numsum(map1_filter.fullMap>0)

            # select points with values less than 0.3
            sum_freq = 0.0
            for fr in range(len(freq_test)):
                sum_freq += float(freq_test[fr])
                # a minimum of 5% (of box size) points are chosen
                if sum_freq/map1_filled > 0.05 and bin_test[fr+1] >= 0.3:
                    t1 = bin_test[fr+1]
                    break
                # if number of points are more than 5% and still have values less than 0.3, a maximum limit of 10% is applied
                if sum_freq/map1_filled > 0.10 or sum_freq > 200000:
                    t1 = bin_test[fr+1]
                    break
            # for the second map
            sum_freq = 0.0
            freq_test = histogram(map2_filter.fullMap,bin_test)[0]
            map2_filled = numsum(map2_filter.fullMap>0)
            for fr in range(len(freq_test)):
                sum_freq += float(freq_test[fr])
                if sum_freq/map2_filled > 0.05 and bin_test[fr+1] >= 0.3:
                    t2 = bin_test[fr+1]
                    break
                if sum_freq/map2_filled > 0.10 or sum_freq > 200000:
                    t2 = bin_test[fr+1]
                    break
            # t1 and t2 are the selected levels based on filtered values and percent of points
            target_points = argwhere((map1_filter.fullMap > 0.0) & (map1_filter.fullMap <= t1))
            probe_points = argwhere((map2_filter.fullMap > 0.0) & (map2_filter.fullMap <= t2))
            if len(target_points) == 0 or len(probe_points) == 0:
                print('Surface detection failed (Mean filter), exiting..')
                return None
            #print len(target_points), len(probe_points), t1, t2
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return None
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return None
        elif Filter == 'Minimum':
            map1_surface = map_target._surface_minimum_filter(float(map_target_threshold1))
            map2_surface = map_probe._surface_minimum_filter(float(map_probe_threshold))
            # select the surface points represented by the mask
            target_points = argwhere(map1_surface == 1)
            probe_points = argwhere(map2_surface == 1)
            if len(target_points) == 0 or len(probe_points) == 0:
                print('Surface detection failed (Minimum filter), exiting..')
                return None
            #print len(target_points), len(probe_points)
            # stop if the number of points are large
            if len(target_points) + len(probe_points) > 250000: return None
            # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return None
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return None
        # surface based on sobel filter on contoured map, high gradient points chosen
        elif Filter == 'Sobel':
            map1_surface = map_target._sobel_filter_contour(float(map_target_threshold1))
            map2_surface = map_probe._sobel_filter_contour(float(map_probe_threshold))

            target_points = argwhere(map1_surface.fullMap > map1_surface.max()/float(2))
            probe_points = argwhere(map2_surface.fullMap > map2_surface.max()/float(2))
            if len(target_points) == 0 or len(probe_points) == 0:
                print('Surface detection failed (Sobel filter), exiting..')
                return None
            #print len(target_points), len(probe_points)
             # check whether the probe points is larger than the probe surface points. if not use the smaller one as probe point
            if len(target_points) < len(probe_points):
                probe_points1 = npcopy(target_points)
                target_points = npcopy(probe_points)
                probe_points = npcopy(probe_points1)
            try:
                from scipy.spatial import cKDTree
                try: kdtree = cKDTree(target_points)
                except RuntimeError: return None
            except ImportError:
                try: kdtree = KDTree(target_points)
                except RuntimeError: return None
        distances = kdtree.query(probe_points)[0]
        #print distances
        #print npmean(distances)
        # by default return mean distance, 1/npmean(distances) gives a similarity score
        if len(distances) == 0: return None
        if not weight:
            if not npmean(distances) <= 0.05: return 1/npmean(distances)
            # becomes inf if mean(dist) is 0. Max score of 20 (will be changed later)
            else: return 1/0.05
            
        x = int(30.0/map_target.apix) # 40A selected as potential distance threshold to calculate weighted score
        if amin(distances) < x/2: distances = distances - amin(distances)
        bins = []
        # to select points that are aligned to target
        i = 0
        while i <= float(x):
            bins.append(i*1.0)
            i += 1
        num_distances = len(distances)
        
        overlap_freq = histogram(distances,bins)[0]
        for fr_i in range(len(overlap_freq)):
            if overlap_freq[fr_i] > amax(overlap_freq)/3.:
                break
            
        total_ext = fr_i
        #might help in accounting for contour difference    
        bins = bins[fr_i:]
        #distancebin = distances < int(x/2.)
        #to check if the aligned surfaces of maps form patches
        if cl:
            ## use this routine to check if the points form a patch
            #points_cl = probe_points[distancebin]
            points_cl = probe_points
            # points_cl represents indices of the smaller map which aligns well with the other map
            # create a kdtree to check whether the points form a patch
            if len(points_cl) == 0: return None,None
            try: kdtree = cKDTree(points_cl)
            except: return None,None
            #cKDtree count_neighbors would work better, but not available in old scipy version
            neighbors_num = 20
            distance_lim = 3.0            
            # query against the same points to check integrity
            neigh = kdtree.query(points_cl,k=neighbors_num,distance_upper_bound=distance_lim)[1]
            ct_neigh = 0
            # for those points where 8 neighbors are not found, len(neigh) is returned as index
            #cl_weight = numsum(numsum(neigh<len(neigh),axis=1) > 15)/float(len(neigh))
            # ratio of 'patch-like' aligned points to total query points : gives the fraction of surface overlap 
            
            cl_weight = numsum(numsum(neigh<len(neigh),axis=1) > 17)/float(len(probe_points))
            # to calculate distances involving these points
            #distances_align = distances[distancebin]
            distances_align = distances
            distances_sel = distances_align[numsum(neigh<len(neigh),axis=1) > 17]
            distances = distances_sel[:]
        
        overlap_freq = histogram(distances,bins)[0]
        total = total_ext #make total_ext=0.0 above for proper contours 
        cumul_freq = 0.0
        enter = 0
        sum_sc = 0.0
        for i in range(len(overlap_freq)):
            w = len(overlap_freq)-(i)
            try:
                cumul_freq += overlap_freq[i]
            except IndexError: pass
            try:
                perc_equiv = float(cumul_freq)/num_distances #/len(distances)
            except ZeroDivisionError:
                print('Distance weighting failed!!. Check surface defined')
                return None, None
            #sum_sc = sum_sc + (npexp(w/2.)*perc_equiv)
            #total += npexp(w/2.)
            sum_sc = sum_sc + ((w)*perc_equiv)
            total += (w)
            enter = 1
        score = float(sum_sc)/total
        if cl:
            if enter == 1:
                if len(distances_sel) == 0.0: return 0.0
                if npmean(distances_sel) == 0.0: return 0.0
                if cl_weight == 0.0: return 0.0
                return score#cl_weight*(1/npmean(distances_sel))
            else: return None, None
        if enter == 1:
            if npmean(distances) <= 0.05: return 1.0
            if npmean(distances) == 0.0: return 1.0
            return score        
        else: return None, None


    def envelope_score(self,map_target, primary_boundary, structure_instance,norm=True):
        """
        
        Calculate the envelope score between a target Map and a Structure Instances.
        
                
        Arguments:
            *map_target*
                Target Map Instance.
            *primary_boundary* 
                Value specified is calculated with primary_boundary of the map object.
            *structure_instance*
                Structure Instance to compare.
        Return:
            Envelope score
            
        """
        
        binMap = map_target.make_bin_map(primary_boundary)
        max_score = float(-2*numsum(binMap.fullMap))
        min_score = float(numsum(binMap.fullMap)-2*numsum(binMap.fullMap+1))
    
        blurrer = StructureBlurrer()
        struct_binMap = blurrer.make_atom_overlay_map1(map_target, structure_instance)
        grid = struct_binMap.get_pos(0.9,1.1)
        for x,y,z in grid:
            g = binMap[z][y][x]
            if g == -1:
                binMap[z][y][x] = 2
            elif g == 0:
                binMap[z][y][x] = -2
        #score=binMap.fullMap.sum()
        score = float(numsum(binMap.fullMap))
        if norm:
            norm_score = float((score-min_score)/(max_score-min_score))
            return norm_score
        else:
            return score

    def envelope_score_map(self,map_target, map_probe,map_target_threshold=0,map_probe_threshold=0,norm=True):
        """
        
        Calculate the envelope score between two Map instance using numoy array. 
        
         Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                EMMap threshold 
                use calcualte_map_threshold to calculate map_target_threshold and map_probe_threshold.                
        Return:
            Envelope score
            
        """
        
        if self.mapComparison(map_target, map_probe):
            if map_target_threshold==0:
                map_target_threshold=self.calculate_map_threshold(map_target)
            if map_probe_threshold==0:
                map_probe_threshold=self.calculate_map_threshold(map_probe)

        binMap = map_target.make_bin_map(map_target_threshold)
        max_score = float(-2*numsum(binMap.fullMap))
        min_score = float(numsum(binMap.fullMap)-2*numsum(binMap.fullMap+1))
        struct_binMap = map_probe.make_bin_map(map_probe_threshold)
        newMap=binMap.fullMap+2*struct_binMap.fullMap
        hist_array=histogram(newMap,4)
        score=2*hist_array[0][0]-(2*(hist_array[0][1]))-(hist_array[0][2])
        #print score, max_score, min_score, numsum(binMap.fullMap)
        if norm:
            norm_score = float((score-min_score))/(max_score-min_score)
            return norm_score
        else:
            return score

        #calculate percent of overlap for two contoured maps
    def _percent_overlap(self,map_target,map_probe,map_target_threshold,map_probe_threshold,flagsize=0):
        """
        
        Calculate the fraction of overlap between two map grids. 
        
         Arguments:
            *map_target, map_probe*
                EMMap instances to compare.
            *map_target_threshold,map_probe_threshold*
                map contour thresholds for map_target and map_probe.                
        Return:
            Percent overlap with respect to smaller grid
            
        """
        if self.mapComparison(map_target,map_probe):
            # contour the first map
            binmap1 = map_target.fullMap > float(map_target_threshold)
            binmap2 = map_probe.fullMap > float(map_probe_threshold)
            # percent calculated on the smaller contoured volume (can be changed)
            minim = len(map_target.fullMap[binmap1])
            if len(map_probe.fullMap[binmap2]) < minim: minim = len(map_probe.fullMap[binmap2])
            maskmap = (binmap1*binmap2) > 0
            if flagsize == 1: return numsum(maskmap), numsum(binmap1), numsum(binmap2)   
            #print numsum(binmap1),numsum(binmap2),numsum(maskmap),minim
            if not minim == 0.0: return float(len(map_target.fullMap[maskmap]))/minim
            else:
                print("Check map contour!!")
                return 0.0
        else:
            print("@@@ Maps could not be matched")
            return -1.0

    def SCCC(self,map_target,resolution_densMap,sigma_map,structure_instance,rigid_body_structure,write=False,c_mode=True):
        """
        
        Calculate Segment based cross-correlation from Pandurangan et al. 2013,J Struct Biol. 2013 Dec 12
        It is a local CCC around a selection of atoms.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
            *structure_instance*
                Structure instance to compare
            *rigid_body_structure*
                Rigid-body Structure instance.
.        Return:
            SCCC score
                
        """
        blurrer = StructureBlurrer()
        scorer = ScoringFunctions()
        outline = ""
        resolution_densMap=float(resolution_densMap)
        whole_fit_map = blurrer.gaussian_blur(structure_instance, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        sim_map = blurrer.gaussian_blur(rigid_body_structure, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        minDens = sim_map.std()
        sim_mask_array = sim_map._get_maskArray(minDens)
        #Apply the mask to em and simulated maps
        mask_emMap=map_target._get_maskMap(sim_mask_array)
        mask_simMap = whole_fit_map._get_maskMap(sim_mask_array)
        #sse_lccf=scorer.CCC(mask_emMap,mask_simMap)
        sse_lccf,ov=scorer.CCC_map(mask_emMap,mask_simMap,cmode=c_mode)
            #return the overall score
        if write==True:
            outline+='SCCC for segment %f\n'%(sse_lccf)
            return outline
        return sse_lccf


    def SCCC_LAP(self,map_target,resolution_densMap,sigma_map,structure_instance,rigid_body_structure,write=False):
        """
        
        Calculate Segment based cross-correlation from Pandurangan et al. 2013,J Struct Biol. 2013 Dec 12
        It is a local CCC around a selection of atoms.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
            *structure_instance*
                Structure instance to compare
            *rigid_body_structure*
                Rigid-body Structure instance.
.        Return:
            SCCC score
                
        """
        blurrer = StructureBlurrer()
        scorer = ScoringFunctions()
        outline = ""
        resolution_densMap=float(resolution_densMap)
        whole_fit_map = blurrer.gaussian_blur(structure_instance, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        sim_map = blurrer.gaussian_blur(rigid_body_structure, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        minDens = sim_map.std()
        sim_mask_array = sim_map._get_maskArray(minDens)
        #Apply the mask to em and simulated maps
        mask_emMap=map_target._get_maskMap(sim_mask_array)
        mask_simMap = whole_fit_map._get_maskMap(sim_mask_array)
        sse_lccf=scorer.laplace_CCC(mask_emMap,mask_simMap)
            #return the overall score
        if write==True:
            outline+='SCCC for segment %f\n'%(sse_lccf)
            return outline
        return sse_lccf


    def SCCC_MI(self,map_target,resolution_densMap,sigma_map,structure_instance,rigid_body_structure,write=False):
        """
        
        Calculate Segment based cross-correlation from Pandurangan et al. 2013,J Struct Biol. 2013 Dec 12
        It is a local CCC around a selection of atoms.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
            *structure_instance*
                Structure instance to compare
            *rigid_body_structure*
                Rigid-body Structure instance.
.        Return:
            SCCC score
                
        """
        blurrer = StructureBlurrer()
        scorer = ScoringFunctions()
        outline = ""
        resolution_densMap=float(resolution_densMap)
        whole_fit_map = blurrer.gaussian_blur(structure_instance, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        sim_map = blurrer.gaussian_blur(rigid_body_structure, resolution_densMap, densMap=map_target, sigma_coeff=sigma_map, normalise=True)
        minDens = sim_map.std()
        sim_mask_array = sim_map._get_maskArray(minDens)
        #Apply the mask to em and simulated maps
        mask_emMap=map_target._get_maskMap(sim_mask_array)
        mask_simMap = whole_fit_map._get_maskMap(sim_mask_array)
        sse_lccf=scorer.MI(mask_emMap,mask_simMap)
            #return the overall score
        if write==True:
            outline+='SCCC for segment %f\n'%(sse_lccf)
            return outline
        return sse_lccf

    def calc_moc(self,indices,map_probe,map_target):
        map_target_mask = map_target.fullMap[indices]
        ##map_target_mask = map_target_mask - float(map_target_mask.sum()/len(map_target_mask))
        map_probe_mask = map_probe.fullMap[indices]
        ##map_probe_mask = map_probe_mask - float(map_probe_mask.sum()/len(map_probe_mask))
        num = numsum(map_target_mask * map_probe_mask)
        den = sqrt(numsum(square(map_target_mask))*numsum(square(map_probe_mask)))
        if den == 0.0: return -1.0
        return num/den

    def SMOC(self,map_target,resolution_densMap,structure_instance,win=11,rigid_body_file=None,sigma_map=0.225,write=False,c_mode=True):
        """
        
        Calculate Local cross correlation (Mander's Overlap)
        It is a local Overlap Coefficient calculated on atoms in sliding residue windows along the chain.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map.
            *structure_instance*
                Model structure instance.
            *win*
                Overlapping Window length to calculate the score
            *rigid_body_file*
                Rigid-body file. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
        Return:
            Dictionary of smoc scores for residues in the chain
        """
        blurrer = StructureBlurrer()
        sim_map = blurrer.gaussian_blur_real_space(structure_instance, resolution_densMap,densMap=map_target,sigma_coeff=sigma_map,normalise=True)
        peak,ave,sigma = sim_map._peak_density()
        #NOTE: filter background
        sim_map.fullMap = sim_map.fullMap*(sim_map.fullMap > peak)
        dict_chain_indices,dict_chain_res, dict_res_dist = blurrer.get_indices(structure_instance,map_target,resolution_densMap,sigma_map)
        #get details of map
        origin = map_target.origin
        apix = map_target.apix
        box_size = map_target.box_size()
        nz,ny,nx = map_target.fullMap.shape
        zg,yg,xg = mgrid[0:nz,0:ny,0:nx]
        indi = list(zip(xg.ravel(), yg.ravel(), zg.ravel()))

        #save rigid body details
        dict_rf_res = {}
        dict_rf_sc = {}
        res_list = []
        rb_list = []
        list_sccc = []
        #save scores for each chain and res
        dict_chain_scores = {}
        #TODO: add multi-chain rigid body parser below
        '''
        r_ct = 0
        if rigid_body_file != None:
            inp = open(rigid_body_file,'r')
            for l in inp:
                if l[0] != '#':
                    score_indices = []
                    lrb = l.split()
                    if len(lrb) == 0: continue
                    r_ct += 1
                    res_list = []
                    rb_pairs = []
                    # get scores for each res and each rigid body
                    for i in range(max((len(lrb)/2)-1,1)):
                        rb_pairs.append([int(lrb[2*i]),int(lrb[2*i+1])])
                        # NOTE: wont work for insertion codes
                        for r in range(int(lrb[2*i]),int(lrb[2*i+1])+1):
                            score_indices.extend(dict_res_indices[r])
                            res_list.append(r)
                    rb_list.append(lrb)
                    dict_rf_res[r_ct] = rb_pairs
                    if len(score_indices) == 0:
                        dict_rf_sc[r_ct] = 0.0#-0.99
                        for res in res_list: dict_res_scores[res] = 0.0#-0.99
                        continue

                    tmplist = score_indices[:]
                    setlist = set(tmplist)
                    score_indices = list(setlist)
                    sc_indices = []
                    for ii in score_indices: sc_indices.append(indi[ii])
                    array_indices = nparray(sc_indices)
                    ind_arrxyz = transpose(array_indices)
                    # get indices for use with map arrays: ([z...],[y...],x...])
                    ind_arrzyx = (ind_arrxyz[2],ind_arrxyz[1],ind_arrxyz[0])
                    sccc = self.calc_moc(ind_arrzyx,sim_map,map_target)
                    dict_rf_sc[r_ct] = sccc
                    #save scores 
                    for res in res_list:
                        dict_res_scores[res] = sccc
                        list_sccc.append(sccc)
            inp.close()
        '''
        #for residues not in rigid bodies: consider pentapeptides
        for ch in dict_chain_indices:
            dict_res_scores = {}
            dict_res_indices = dict_chain_indices[ch]    
            for res in dict_res_indices:
                if not dict_res_scores.has_key(res):
                    indices = dict_res_indices[res][:]
                    #consider residues on both sides. NOTE: wont work for insertion codes!
                    #need to rewite res numbers to avoid insertion codes
                    for ii in range(1,int(round((win+1)/2))):
                        try:
                            #get prev residue indices
                            indices.extend(dict_res_indices[dict_chain_res[ch][dict_chain_res[ch].index(res)-ii]])
                        except: pass
                    for ii in range(1,int(round((win+1)/2))):
                        try:
                            indices.extend(dict_res_indices[dict_chain_res[ch][dict_chain_res[ch].index(res)+ii]])
                        except: pass

                    tmplist = indices[:]
                    setlist = set(tmplist)
                    indices = list(setlist)
                    sc_indices = []
                    for ii in indices: sc_indices.append(indi[ii])
                    if len(indices) < 10:
                        try: 
                            dict_res_scores[res] = dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)-1]]
                            try: dict_res_scores[res] = (dict_res_scores[res]+dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)+1]])/2.0
                            except (IndexError,KeyError): pass
                        except (IndexError,KeyError): 
                            try: dict_res_scores[res] = dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)+1]]
                            except (IndexError,KeyError): dict_res_scores[res] = 0.0
                        continue
                    array_indices = nparray(sc_indices)
                    ind_arrxyz = transpose(array_indices)
                    ind_arrzyx = (ind_arrxyz[2],ind_arrxyz[1],ind_arrxyz[0])
                    sccc = self.calc_moc(ind_arrzyx,sim_map,map_target)
                    dict_res_scores[res] = sccc
                    if sccc == -1.0:
                        try:
                            dict_res_scores[res] = dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)-1]]
                            try: dict_res_scores[res] = (dict_res_scores[res]+dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)+1]])/2.0
                            except (IndexError,KeyError): pass
                        except IndexError:
                            try: dict_res_scores[res] = dict_res_scores[dict_chain_res[ch][dict_chain_res[ch].index(res)+1]]
                            except (IndexError,KeyError): dict_res_scores[res] = 0.0
                        continue
                    list_sccc.append(sccc)
            dict_chain_scores[ch] = dict_res_scores
        return dict_chain_scores, dict_chain_res


    def _SMOC1(self,map_target,resolution_densMap,structure_instance,win=11,rigid_body_file=None,sigma_map=0.225,write=False):
        """
        
        Calculate Local cross correlation (Mander's Overlap)
        It is a local Overlap Coefficient calculated on atoms in sliding residue windows along the chain.  
                
        Arguments:
        
            *map_target*
                Target Map Instance.
            *resolution_densMap*
                Parameter need for Structure Blurrer.
                Resolution of the target map.
            *structure_instance*
                Model structure instance.
            *win*
                Overlapping Window length to calculate the score
            *rigid_body_file*
                Rigid-body file. 
            *sigma_map*
                Parameter need for Structure Blurrer.
                The sigma value (multiplied by the resolution) that controls the width of the Gaussian. 
                Default values is 0.356.
                
                Other values used :
                
                    0.187R corresponding with the Gaussian width of the Fourier transform falling to half the maximum at 1/resolution, as used in Situs (Wriggers et al, 1999);
                    
                    0.356R corresponding to the Gaussian width at 1/e maximum height equaling the resolution, the default in Chimera (Petterson et al, 2004);
                    
                    0.425R the fullwidth half maximum being equal to the resolution, as used by FlexEM (Topf et al, 2008);
                                
                    0.5R the distance between the two inflection points being the same length as the resolution, an option in Chimera (Petterson et al, 2004);
                                
                    1R where the sigma value simply equal to the resolution, as used by NMFF (Tama et al, 2004).
        Return:
            Dictionary of smoc scores for residues in the chain
        """
         
        blurrer = StructureBlurrer()
        sim_map = blurrer.gaussian_blur_real_space(structure_instance, resolution_densMap,densMap=map_target,sigma_coeff=sigma_map,normalise=True)
        peak,ave,sigma = sim_map._peak_density()
        #NOTE: filter background
        sim_map.fullMap = sim_map.fullMap*(sim_map.fullMap > peak)
        dict_res_indices,dict_res_dist = blurrer.get_indices(structure_instance,map_target,resolution_densMap)
        #get details of map
        origin = map_target.origin
        apix = map_target.apix
        box_size = map_target.box_size()
        nz,ny,nx = map_target.fullMap.shape
        zg,yg,xg = mgrid[0:nz,0:ny,0:nx]
        indi = list(zip(xg.ravel(), yg.ravel(), zg.ravel()))

        #save rigid body details
        dict_rf_res = {}
        dict_rf_sc = {}
        res_list = []
        rb_list = []
        list_sccc = []
        #save scores for each res
        dict_res_scores = {}
        r_ct = 0
        if rigid_body_file != None:
            inp = open(rigid_body_file,'r')
            for l in inp:
                if l[0] != '#':
                    score_indices = []
                    lrb = l.split()
                    if len(lrb) == 0: continue
                    r_ct += 1
                    res_list = []
                    rb_pairs = []
                    # get scores for each res and each rigid body
                    for i in range(max((len(lrb)/2)-1,1)):
                        rb_pairs.append([int(lrb[2*i]),int(lrb[2*i+1])])
                        # NOTE: wont work for insertion codes
                        for r in range(int(lrb[2*i]),int(lrb[2*i+1])+1):
                            score_indices.extend(dict_res_indices[r])
                            res_list.append(r)
                    rb_list.append(lrb)
                    dict_rf_res[r_ct] = rb_pairs
                    if len(score_indices) == 0:
                        dict_rf_sc[r_ct] = 0.0#-0.99
                        for res in res_list: dict_res_scores[res] = 0.0#-0.99
                        continue

                    tmplist = score_indices[:]
                    setlist = set(tmplist)
                    score_indices = list(setlist)
                    sc_indices = []
                    for ii in score_indices: sc_indices.append(indi[ii])
                    array_indices = nparray(sc_indices)
                    ind_arrxyz = transpose(array_indices)
                    # get indices for use with map arrays: ([z...],[y...],x...])
                    ind_arrzyx = (ind_arrxyz[2],ind_arrxyz[1],ind_arrxyz[0])
                    sccc = self.calc_moc(ind_arrzyx,sim_map,map_target)
                    dict_rf_sc[r_ct] = sccc
                    #save scores 
                    for res in res_list:
                        dict_res_scores[res] = sccc
                        list_sccc.append(sccc)
            inp.close()
        #for residues not in rigid bodies: consider pentapeptides    
        for res in dict_res_indices:
            if not dict_res_scores.has_key(res):
                indices = dict_res_indices[res][:]
                #consider residues on both sides. NOTE: wont work for insertion codes!
                #need to rewite res numbers to avoid insertion codes
                for ii in range(1,int(round((win+1)/2))):
                    try:
                        indices.extend(dict_res_indices[res-ii])
                    except: pass
                for ii in range(1,int(round((win+1)/2))):
                    try:
                        indices.extend(dict_res_indices[res+ii])
                    except: pass

                tmplist = indices[:]
                setlist = set(tmplist)
                indices = list(setlist)
                sc_indices = []
                for ii in indices: sc_indices.append(indi[ii])
                if len(indices) == 0:
                    dict_res_scores[res] = 0.0#-0.99
                    continue
                array_indices = nparray(sc_indices)
                ind_arrxyz = transpose(array_indices)
                ind_arrzyx = (ind_arrxyz[2],ind_arrxyz[1],ind_arrxyz[0])
                sccc = self.calc_moc(ind_arrzyx,sim_map,map_target)
                dict_res_scores[res] = sccc
                list_sccc.append(sccc)
        return dict_res_scores

    def _get_shell(self,dist1,maxlevel,step):
        # indices between upper and lower shell bound
        fshells1 = ((dist1 < min(maxlevel,x+step)) & (dist1 >= x))

        # match power spectra for two maps
    def _amplitude_match(self,map_1,map_2,shellmin,shellmax,step=0.005,c1=0,c2=0,reso=None,lpfiltb=False,lpfilta=False,ref=False):
        # fourier transform: use pyfftw if available
        pyfftw_flag = 1
        try:
            import pyfftw
        except ImportError: pyfftw_flag = 0
        try:
            if pyfftw_flag == 0: raise ImportError
            inputa1 = pyfftw.n_byte_align_empty(map_1.fullMap.shape, 16, 'complex128')
            outputa1 = pyfftw.n_byte_align_empty(map_1.fullMap.shape, 16, 'complex128')
            # fft planning, set planning_timelimit or flags to make it faster
            fft = pyfftw.FFTW(inputa1,outputa1,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa1[:,:,:] = map_1.fullMap[:,:,:]
            fft()
            ft1 = Map(fftshift(outputa1), map_1.origin, map_1.apix, map_1.filename, map_1.header[:])
        except:
            # use numpy fft instead
            ft1 = map_1.fourier_transform()
        try:
            if pyfftw_flag == 0: raise ImportError
            inputa2 = pyfftw.n_byte_align_empty(map_2.fullMap.shape, 16, 'complex128')
            outputa2 = pyfftw.n_byte_align_empty(map_2.fullMap.shape, 16, 'complex128')
            fft = pyfftw.FFTW(inputa2,outputa2,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa2[:,:,:] = map_2.fullMap[:,:,:]
            fft()
            ft2 = Map(fftshift(outputa2), map_2.origin, map_2.apix, map_2.filename, map_2.header[:])
        except:
            ft2 = map_2.fourier_transform()
        #low pass filter before scaling
        if reso != None:
            cutoff1 = map_1.apix/float(reso)
            cutoff2 = map_2.apix/float(reso)
            if lpfiltb and not lpfilta:
                ft1._tanh_lowpass(cutoff1,fall=0.2,ftmap=True)
                ft2._tanh_lowpass(cutoff2,fall=0.2,ftmap=True)
        # max dimension
        size1 = max(map_1.x_size(),map_1.y_size(),map_1.z_size())
        #shell values correspond to freq: 0-0.5 (nyquist)
        dist1 = map_1._make_fourier_shell(1)/map_1.apix
        size2 = max(map_2.x_size(),map_2.y_size(),map_2.z_size())
        #shell values correspond to freq: 0-0.5 (nyquist)
        dist2 = map_2._make_fourier_shell(1)/map_2.apix
        #SCALING
        # storing for plots
        ft1_avg = []
        ft2_avg = []
        ft1_avg_new = []
        lfreq = []
        # select max spatial frequency to iterate to. low resolution map
        maxlevel = 0.5/max(map_1.apix,map_2.apix)
        # loop over freq shells, shellwidth=0.005
        #for x in arange(0,maxlevel+step,step):
        nc = 0
        x = 0.0
        highlevel = x+step
        while (x<maxlevel):
            #print x,highlevel, maxlevel
            # indices between upper and lower shell bound
            fshells1 = ((dist1 < min(maxlevel,highlevel)) & (dist1 >= x))
            # radial average
            shellvec1 = ft1.fullMap[fshells1]
            # indices between upper and lower shell bound
            fshells2 = ((dist2 < min(maxlevel,highlevel)) & (dist2 >= x))
            # radial average
            shellvec2 = ft2.fullMap[fshells2]
            #if len(shellvec1) == 0 or len(shellvec2) == 0: continue    
            abs1 = abs(shellvec1)
            abs2 = abs(shellvec2)
            #print nonzero(abs1)
            #print nonzero(abs2)
            ns1 = len(nonzero(abs1)[0]) #or count_nonzero
            ns2 = len(nonzero(abs2)[0]) #or count_nonzero
            if ns1 < 10 or ns2 < 10:
                nc += 1
                highlevel = min(maxlevel,x+(nc+1)*step) 
                x = max(0.0,x-nc*step)
                #print ns1, ns2
                continue
            else: nc = 0
                    
            mft1 = npmean(abs1)#npmean(sqrt(shellvec1.real**2+shellvec1.imag**2))
            mft2 = npmean(abs2)#npmean(sqrt(shellvec2.real**2+shellvec2.imag**2))#npmean(abs(ft2.fullMap[fshells2]))
            if mft1 == 0.0 and mft2 == 0.0:
                continue
            # sq of radial avg amplitude
            ft1_avg.append(np_log10(npmean(square(abs1))))
            ft2_avg.append(np_log10(npmean(square(abs2))))
            

            # scale to amplitudes of the ref map
            if ref:
                if mft1 == 0.0: continue
                ft1.fullMap[fshells1] = shellvec1*(mft2/mft1)
            else:
                # replace with avg amplitudes for the two maps
                ft1.fullMap[fshells1] = shellvec1*(mft2+mft1)/(2*mft1)
                ft2.fullMap[fshells2] = shellvec2*(mft2+mft1)/(2*mft2)

            # new radial average (to check)
            mft1 = npmean(abs(ft1.fullMap[fshells1]))#numsum(absolute(ft1.fullMap[fshells1]))/len(shellvec1)
            ft1_avg_new.append(np_log10(npmean(square(abs(ft1.fullMap[fshells1])))))
            lfreq.append(highlevel)

            sampling_frq = highlevel

            cutoff_freq = min((1.0/reso) + 0.25,maxlevel) # 0.25 added to reso based cutoff
            #print 'freq cutoff', (1.0/reso)+0.25, maxlevel
            
            # scale the rest and break after relevant frequencies
            if sampling_frq > cutoff_freq:
                fshells1 = (dist1 >= highlevel)
                shellvec1 = ft1.fullMap[fshells1]
                mft1 = npmean(abs(shellvec1))
                fshells2 = (dist2 >= highlevel)
                shellvec2 = ft2.fullMap[fshells2]
                mft2 = npmean(abs(shellvec2))
                if mft1 == 0.0 and mft2 == 0.0:
                    break
                ft1_avg.append(np_log10(npmean(square(abs(shellvec1)))))
                ft2_avg.append(np_log10(npmean(square(abs(shellvec2)))))

                if ref:
                    if mft1 == 0.0: break
                    ft1.fullMap[fshells1] = shellvec1*(mft2/mft1)
                else:
                    ft1.fullMap[fshells1] = shellvec1*(mft2+mft1)/(2*mft1)
                    ft2.fullMap[fshells2] = shellvec2*(mft2+mft1)/(2*mft2)

                mft1 = npmean(abs(ft1.fullMap[fshells1])) #after scaling
                ft1_avg_new.append(np_log10(npmean(square(abs(ft1.fullMap[fshells1]))))) #after scaling
                lfreq.append((highlevel+step/2))
                break
            x = highlevel
            highlevel = x+step
        # low pass filter after?
        #low pass filter before scaling
        if reso != None:
            if lpfilta and not lpfiltb:
                ft1._tanh_lowpass(cutoff1,fall=0.2,ftmap=True)
                ft2._tanh_lowpass(cutoff2,fall=0.2,ftmap=True)

        # ifft
        try:
            if pyfftw_flag == 0: raise ImportError
            ifft = pyfftw.FFTW(inputa1,outputa1,direction='FFTW_BACKWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa1[:,:,:] = ifftshift(ft1.fullMap)[:,:,:]
            ifft()
            map1_filt = Map(outputa1.real.astype('float'), map_1.origin, map_1.apix, map_1.filename, map_1.header[:])
        except:
            # use numpy ifft instead
            map1_filt = map_1.copy()
            map1_filt.fullMap = real(ifftn(ifftshift(ft1.fullMap)))
        try:
            if pyfftw_flag == 0: raise ImportError
            ifft = pyfftw.FFTW(inputa2,outputa2,direction='FFTW_BACKWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa2[:,:,:] = ifftshift(ft2.fullMap)[:,:,:]
            ifft()
            map2_filt = Map(outputa2.real.astype('float'), map_2.origin, map_2.apix, map_2.filename, map_2.header[:])
        except:
            map2_filt = map_2.copy()
            map2_filt.fullMap = real(ifftn(ifftshift(ft2.fullMap)))

        
        try:
            # to check frequency plots
            #print lfreq
            #print ft1_avg
            #print ft2_avg
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from matplotlib import pylab
            try: plt.style.use('ggplot')
            except AttributeError: pass
            plt.rcParams.update({'font.size': 18})
            plt.rcParams.update({'legend.fontsize': 18})
            plt.plot(lfreq,ft1_avg,'r--',label='map1')
            plt.plot(lfreq,ft2_avg,'bs',label='map2')
            plt.plot(lfreq,ft1_avg_new,'g^',label='scaled')
            #plt.show()
            leg = plt.legend(loc='upper right')
            for legobj in leg.legendHandles:
                legobj.set_linewidth(2.0)
            pylab.savefig("spectra.png")
            plt.close()
        except: pass
        
        return map1_filt.fullMap,map2_filt.fullMap

        # FSC for two maps
    def _fsc(self,map_1,map_2,shellmin,shellmax,step=0.005,c1=0,c2=0,reso=None):
        # fourier transform: use pyfftw if available
        pyfftw_flag = 1
        try:
            import pyfftw
        except ImportError: pyfftw_flag = 0
        try:
            if pyfftw_flag == 0: raise ImportError
            inputa1 = pyfftw.n_byte_align_empty(map_1.fullMap.shape, 16, 'complex128')
            outputa1 = pyfftw.n_byte_align_empty(map_1.fullMap.shape, 16, 'complex128')
            # fft planning, set planning_timelimit or flags to make it faster
            fft = pyfftw.FFTW(inputa1,outputa1,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa1[:,:,:] = map_1.fullMap[:,:,:]
            fft()
            ft1 = Map(fftshift(outputa1), map_1.origin, map_1.apix, map_1.filename, map_1.header[:])
        except:
            # use numpy fft instead
            ft1 = map_1.fourier_transform()
        try:
            if pyfftw_flag == 0: raise ImportError
            inputa2 = pyfftw.n_byte_align_empty(map_2.fullMap.shape, 16, 'complex128')
            outputa2 = pyfftw.n_byte_align_empty(map_2.fullMap.shape, 16, 'complex128')
            fft = pyfftw.FFTW(inputa2,outputa2,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa2[:,:,:] = map_2.fullMap[:,:,:]
            fft()
            ft2 = Map(fftshift(outputa2), map_2.origin, map_2.apix, map_2.filename, map_2.header[:])
        except:
            ft2 = map_2.fourier_transform()
        #low pass filter before scaling
        if reso != None:
            cutoff1 = map_1.apix/float(reso)
            cutoff2 = map_2.apix/float(reso)
        # max dimension
        size1 = max(map_1.x_size(),map_1.y_size(),map_1.z_size())
        #shell values correspond to freq: 0-0.5 (nyquist)
        #and convert to abs frequencies
        dist1 = map_1._make_fourier_shell(1)/map_1.apix
        size2 = max(map_2.x_size(),map_2.y_size(),map_2.z_size())
        #SCALING
        # storing for plots
        lfreq = []
        # select max spatial frequency to iterate to. low resolution map
        maxlevel = 0.5/max(map_1.apix,map_2.apix)
        # loop over freq shells, shellwidth=0.005
        #for x in arange(0,maxlevel+step,step):
        nc = 0
        x = 0.0
        listC = []
        highlevel = x+step
        while (x<maxlevel):
            #print x,highlevel, maxlevel
            # indices between upper and lower shell bound
            C1 = 0.0
            C2 = 0.0
            C3 = 0.0
            fshells = argwhere((dist1 < min(maxlevel,highlevel)) & (dist1 >= x))        
            # shell values
            shellvec1 = ft1.fullMap[transpose(fshells)]
            # shell values
            shellvec2 = ft2.fullMap[transpose(fshells)]
            #if len(shellvec1) == 0 or len(shellvec2) == 0: continue    
            abs1 = abs(shellvec1)
            abs2 = abs(shellvec2)
            #print nonzero(abs1)
            #print nonzero(abs2)
            ns1 = len(nonzero(abs1)[0]) #or count_nonzero
            ns2 = len(nonzero(abs2)[0]) #or count_nonzero
            if ns1 < 10 or ns2 < 10:
                nc += 1
                highlevel = min(maxlevel,x+(nc+1)*step) 
                x = max(0.0,x-nc*step)
                #print ns1, ns2
                continue
            else: nc = 0
            for v in fshells:
                if v[2] > 0 or (v[0] >= 0 and (v[1] >= 0 or v[0] != 0)):
                    C1 += ft1.fullMap[v[0]][v[1]][v[2]]*conjugate(ft2.fullMap[v[0]][v[1]][v[2]])
                    C2 += ft1.fullMap[v[0]][v[1]][v[2]]*conjugate(ft1.fullMap[v[0]][v[1]][v[2]])
                    C3 += ft2.fullMap[v[0]][v[1]][v[2]]*conjugate(ft2.fullMap[v[0]][v[1]][v[2]])
                    
            listC.append(abs(C1)/sqrt(abs(C2)*abs(C3)))
            print(abs(C1)/sqrt(abs(C2)*abs(C3)), (x+highlevel)/2.)
            lfreq.append(highlevel)

            sampling_frq = highlevel

            cutoff_freq = min((1.0/reso) + 0.25,maxlevel) # 0.1 added to reso based cutoff
            #print 'freq cutoff', (1.0/reso), sampling_frq/map_1.apix
            
            # scale the rest and break after relevant frequencies
            if sampling_frq > cutoff_freq:
                fshells = argwhere(dist1 >= highlevel)
                for v in fshells:
                    if v[2] > 0 or (v[0] >= 0 and (v[1] >= 0 or v[0] != 0)):
                        C1 += ft1.fullMap[v[0]][v[1]][v[2]]*conjugate(ft2.fullMap[v[0]][v[1]][v[2]])
                        C2 += ft1.fullMap[v[0]][v[1]][v[2]]*conjugate(ft1.fullMap[v[0]][v[1]][v[2]])
                        C3 += ft2.fullMap[v[0]][v[1]][v[2]]*conjugate(ft2.fullMap[v[0]][v[1]][v[2]])
                    
                listC.append(abs(C1)/sqrt(abs(C2)*abs(C3)))
                print(abs(C1)/sqrt(abs(C2)*abs(C3)), (x+highlevel)/2.)
            
                lfreq.append((highlevel+step/2))
                break
            x = highlevel
            highlevel = x+step
        
        # to check frequency plots
        import matplotlib
        #matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import pylab
        fig = plt.plot(lfreq,listC,'g^')
        plt.show()
        pylab.savefig("test.png")
        plt.close()
        
        return
        #Guess not requited here. Check and remove it.
    def get_clash_map(self,emmap, apix):
        template_grid = emmap._make_clash_map(apix)
        return template_grid

    def get_sm_score(self, struct, ncomp, template_grid, cvol, apix):
        overlay_maplist = []
        overlay_maplist = self.get_overlay_comp_maplist(struct,template_grid)
        nc = list(range(ncomp))
        cpair = list(itertools.combinations(nc,2))
        score = 0.0
        n_overlap_voxel = 0
        overlap_volume = 0.0
        for i in cpair:
                n_overlap_voxel = (overlay_maplist[i[0]].fullMap * overlay_maplist[i[1]].fullMap).sum()
                #overlap_volume = ((n_overlap_voxel*2)*apix)**3
                overlap_volume = ((apix**3)*n_overlap_voxel) * 2
                clash_percent = (float(overlap_volume / (cvol[i[0]]+cvol[i[1]])))
                score = score + clash_percent
        return -(score)

    def get_overlay_comp_maplist(self, struct,template_grid):
        blurrer = StructureBlurrer()
        overlay_maplist = []
        #ssplit = struct.structList
        ssplit = struct.split_into_chains()
        #split_into_chains()
        for x in ssplit:
                #print 'Chain:'
                #CHANGE HERE FOR SHAPE SCORE BASED ON OVERLAP SCORE OR GRID SCORE
                overlay_maplist.append(blurrer.make_atom_overlay_map1(template_grid, x))
                #print 'chain ids from overlay_maplist ', x
                #overlay_maplist.append(blurrer.get_shapeGrid(template_grid, x))

        #print 'Done overlay_comp_maplist'
        #exit(0)
        return overlay_maplist
