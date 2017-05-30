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

from numpy import  array, int32, float32, zeros, real, argwhere, diag, histogram, dot, matrix, amin, arange,\
                   indices, ravel, all as all_points, delete, transpose, searchsorted, newaxis, where, meshgrid,\
                   ma, sum as numsum,median,sqrt as srt, digitize, nonzero, floor, ceil, amax, mean as npmean,\
                   std as npstd, square as npsquare, tanh as np_tanh,set_printoptions
from random import randrange
from scipy.ndimage.interpolation import  shift, affine_transform, zoom, map_coordinates, spline_filter
from scipy.ndimage import laplace, uniform_filter, generic_filter,minimum_filter
from scipy.ndimage.filters import sobel
from scipy.fftpack import fftn, ifftn, fftshift,fftfreq, ifftshift
from scipy.signal import resample
from scipy.spatial import KDTree
from scipy import weave
from scipy.weave import converters
from scipy.ndimage.morphology import binary_opening
from scipy.ndimage import measurements
import sys, datetime
import struct as binary
import TEMPy.Vector  as Vector
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom
import math


class Map:
    """ 
    
    A class representing all information from a density map file. 
    NOTE: Currently it can only read the CCP4/MRC  format.
    
    """

    def __init__(self, fullMap, origin, apix, filename, header=[]):
        """
        
        Read a map and its parameters in to Map class instance.
        
        *filename*
            name of map file.
        *origin*    
            origin co-ordinates of the map (x_origin, y_origin, z_origin).
        *apix*
            grid spacing of map.
        *filename*
            filename of the Map instance
            
            NOTE: The *filename* 'build++copy' is reserved for copying of other Map class instances."""        
        self.header = header
        self.origin = origin
        self.apix = apix
        self.filename = filename
        self.fullMap = fullMap

    def __repr__(self):
        box_size = list(self.box_size())
        box_size.reverse()
        string1 = 'Obtained from ' + self.filename + '\n'
        string2 = 'Origin: '+ str(self.origin) + '\n'
        string3 = 'Box size (x,y,z): ' + str(box_size) + '\n'
        string4 = 'Grid spacing: ' + str(self.apix) + '\n'
        string5 = 'Min, max, mean, std: %.3f, %.3f, %.3f, %.3f \n' %(self.min(), self.max(), self.mean(), self.std())
        return string1 + string2 + string3 + string4 +string5

    def x_origin(self):
        """
         Return:
             x-coordinate of the origin.
        
        """
        return self.origin[0]

    def y_origin(self):
        """
        
        Return:
            y-coordinate of the origin.
        
        """
        return self.origin[1]

    def z_origin(self):
        """
        
        Return:
            z-coordinate of the origin.
        
        """
        return self.origin[2]
    
    def copy(self):
        """
        
        Return:
            copy of the Map.
        
        """
        copy = Map(self.fullMap.copy(), self.origin[:], self.apix, self.filename, self.header[:])
        return copy

    def getMap(self):
        """
        
        Return:
            array containing the map density data.
        
        """
        return self.fullMap

    def box_size(self):
        """
        
        Return:
            size of the map array, in ZYX format.
        
        """
        return self.fullMap.shape

    def x_size(self):
        """
        
        Return:
            x size of the map array in x direction.
        
        """
        return self.fullMap.shape[2]
    
    def y_size(self):
        """
        
        Return:
            y size of the map array in y direction.
        
        """
        return self.fullMap.shape[1]
    
    def z_size(self):
        """
        
        Return:
            z size of the map array in z direction.
        
        """
        return self.fullMap.shape[0]
    def map_size(self):
        """
        
        Return:
            size of the array fullMap.
        
        """
        return self.fullMap.size

    def __getitem__(self, index):
        """
        
        Allows direct indexing of the map array from the Map instance.
        ie. map[2][2][1] instead of map.fullMap[2][2][1]
        
        """
        return self.fullMap[index]

    # -- Map modification methods. All of them return a new Map instance. -- #
    def _shift_density(self, offset):
        self.fullMap = self.fullMap + float(offset)

    def scale_map(self, scaling):
        """
        Scaling Map by scaling factor

        Return:
            new Map instance         
         """
        sc = 1./scaling
        c_sh = self.pixel_centre()*(1-sc)
        newMap = self.copy()
        newMap.fullMap = affine_transform(self.fullMap, diag([sc,sc,sc]), offset=[c_sh.z, c_sh.y, c_sh.x])
        return newMap

    def _crop_box(self,c,f):
        """
        Crop a map in place based on a threshold
        Arguments:
                *c*
                        map threshold
                *f*
                        factor to relax threshold
        Return:
                
        """
        minval = float(c) - (float(f)*self.std())
        
        axis_diff = []
        for i in range(3):
            ct1 = 0
            try:
                while (self.fullMap[0] < minval).all(): #numsum(self.fullMap[0] < minval)/self.fullMap[0].size > 0.99
                    self.fullMap = delete(self.fullMap,0,0)
                    ct1 += 1
            except (IndexError, ValueError) as exc: pass
            axis_diff.append(ct1)
            ct2 = 0
            try:
                while (self.fullMap[self.fullMap.shape[0]-1] < minval).all():
                    self.fullMap = delete(self.fullMap,-1,0)
                    ct2 += 1
            except (IndexError, ValueError) as exc: pass
            self.fullMap = transpose(self.fullMap,(2,0,1))
        ox = self.origin[0] + axis_diff[1]*self.apix
        oy = self.origin[1] + axis_diff[2]*self.apix
        oz = self.origin[2] + axis_diff[0]*self.apix
        self.origin = (ox,oy,oz)
        '''
        #ignore background slices based on threshold 
        zi,yi,xi = nonzero(self.fullMap > minval)
        zimin,zimax = zi.min(),zi.max()
        yimin,yimax = yi.min(),yi.max()
        ximin,ximax = xi.min(),xi.max()
        self.fullMap = self.fullMap[zimin:zimax+1,yimin:yimax+1,ximin:ximax+1]
        ox = self.origin[0] + ximin*self.apix
        oy = self.origin[1] + yimin*self.apix
        oz = self.origin[2] + zimin*self.apix
        self.origin = (ox,oy,oz)
        '''
        
    def _alignment_box(self,map2,s):
        m1_shape = self.fullMap.shape
        m2_shape = map2.fullMap.shape
        (ox,oy,oz) = (self.origin[0],self.origin[1],self.origin[2])
        (o1x,o1y,o1z) = (map2.origin[0],map2.origin[1],map2.origin[2])
        offset = (o1x-ox,o1y-oy,o1z-oz)
        (m1x,m1y,m1z) = (ox+self.fullMap.shape[2]*self.apix,oy+self.fullMap.shape[1]*self.apix,oz+self.fullMap.shape[0]*self.apix)
        (m2x,m2y,m2z) = (o1x+map2.fullMap.shape[2]*map2.apix,o1y+map2.fullMap.shape[1]*map2.apix,o1z+map2.fullMap.shape[0]*map2.apix)
        (nx,ny,nz) = (o1x, o1y, o1z)
        if offset[0] > 0:
                nx = ox
        if offset[1] > 0:
                ny = oy
        if offset[2] > 0:
                nz = oz
	
        (lz,ly,lx) = ((m2z-nz)/float(s),(m2y-ny)/float(s),(m2x-nx)/float(s))
        if m2x < m1x: lx = (m1x-nx)/float(s)
        if m2y < m1y: ly = (m1y-ny)/float(s)
        if m2z < m1z: lz = (m1z-nz)/float(s)
        gridshape = (int(lz),int(ly),int(lx))
        new_origin = (nx,ny,nz)
        return gridshape, new_origin

    def _interpolate_to_grid(self,grid,s,ori):
        new_map = self.copy()
        (ox,oy,oz) = (self.origin[0],self.origin[1],self.origin[2])
        (o1x,o1y,o1z) = (float(ori[0]),float(ori[1]),float(ori[2]))
        scale = float(s)/self.apix
        offset = (o1x-ox,o1y-oy,o1z-oz)

        gridshape = grid
        new_map.origin = (o1x, o1y, o1z)
        grid_indices = indices(gridshape)
        z_ind = grid_indices[0]
        z_ind.ravel()
        y_ind = grid_indices[1]
        y_ind.ravel()
        x_ind = grid_indices[2]
        x_ind.ravel()
        z_ind = ((offset[2])/self.apix)+scale*z_ind
        y_ind = ((offset[1])/self.apix)+scale*y_ind
        x_ind = ((offset[0])/self.apix)+scale*x_ind
        new_array = map_coordinates(self.fullMap,[z_ind,y_ind,x_ind],cval=self.min())
        new_map.fullMap = new_array.reshape(gridshape)
        new_map.origin = (o1x, o1y, o1z)
        new_map.apix = float(s)
        return new_map
    '''
    def _grid_footprint(self):
        a = zeros((3,3,3))
        a[1,1,1] = 1
        a[0,1,1] = 1
        a[1,0,1] = 1
        a[1,1,0] = 1
        a[2,1,1] = 1
        a[1,2,1] = 1
        a[1,1,2] = 1
        return a
    '''
    #for maps with different grid spacing along the 3 axes, use this function to downsample
    def _downsample_apix(self,spacing):
        #downsample
        apix_ratio = (round(self.header[10]/self.header[7],2)/spacing,round(self.header[11]/self.header[8],2)/spacing,round(self.header[12]/self.header[9],2)/spacing)
        grid_shape = int(round(self.z_size()*apix_ratio[2])),int(round(self.y_size()*apix_ratio[1])),int(round(self.x_size()*apix_ratio[0]))
        try: newmap = self._interpolate_to_grid1(grid_shape,spacing,self.origin)
        except: newmap = self._interpolate_to_grid(grid_shape,spacing,self.origin)
        return newmap
  
    # find background peak
    def _peak_density(self):
        """
        
        Find background peak and sigma (for values beyond the peak)

        Return:
            peak, average and sigma (beyond peak)         
        
        """        
        freq,bins = histogram(self.fullMap,1000)
        ind = nonzero(freq==amax(freq))[0]
        peak = None
        ave = npmean(self.fullMap)
        sigma = npstd(self.fullMap)
        
        for i in ind:
            val = (bins[i]+bins[i+1])/2.
            if val < float(ave)+float(sigma):
                peak = val
        if peak == None:
            peak = ave
        sigma1 = None
        if peak != None:
            mask_array = self.fullMap[self.fullMap > peak]
            sigma1 = srt(npmean(npsquare(mask_array-peak)))
        
        return peak,ave,sigma1
        

    def _sobel_surface_mask(self,c):
      """
        
      Apply sobel filter on binned density maps 
 
      Return:
            new Map instance         
        
      """

      newmap = self.copy()
      binmap = newmap.fullMap > float(c)
      sx = sobel(binmap,0,mode='constant')
      sy = sobel(binmap,1,mode='constant')
      sz = sobel(binmap,2,mode='constant')
      newmap.fullMap = srt(sx*sx+sy*sy+sz*sz)
      newmap.fullMap = binmap*newmap.fullMap
      return newmap

    def _sobel_filter_contour(self,c):
        """
        
       Apply sobel filter on density maps above contour

        Return:
            new Map instance         
        
        """

        newmap = self.copy()
        
        binmap = newmap.fullMap > c
        # use the line below to apply the filter on a binary contoured map
        newmap.fullMap = binmap*newmap.fullMap
        # apply kernel convolution along each axis
        sx = sobel(newmap.fullMap,0,mode='constant')
        sy = sobel(newmap.fullMap,1,mode='constant')
        sz = sobel(newmap.fullMap,2,mode='constant')
        # calculate gradient magnitude
        newmap.fullMap = srt(sx*sx+sy*sy+sz*sz)
        # to select points only within the contour, use the line below
        ##newmap.fullMap = binmap*newmap.fullMap
        return newmap
        
    def _sobel_filter_map_all(self):
        """
        
       Apply sobel filter on density maps

        Return:
            new Map instance         
        
        """

        newmap = self.copy()
        #binmap = newmap.fullMap > float(c)
        #print binmap
        #newmap.fullMap = binmap*newmap.fullMap
        #print newmap.fullMap
        #newmap.write_to_MRC_file('msobelpre.mrc')
        sx = sobel(newmap.fullMap,0,mode='constant')
        sy = sobel(newmap.fullMap,1,mode='constant')
        sz = sobel(newmap.fullMap,2,mode='constant')
        newmap.fullMap = srt(sx*sx+sy*sy+sz*sz)
        #newmap.write_to_MRC_file('msobelpost.mrc')
        #newmap.fullMap = binmap*newmap.fullMap
        #newmap.write_to_MRC_file('msobelpostBIN.mrc')
        return newmap
 
    def _laplace_filtered_contour(self,c):
        """
        
       Apply Laplacian filter on density maps above contour

        Return:
            new Map instance         
        
        """
        newmap = self.copy()
        binmap = newmap.fullMap > float(c)
        newmap.fullMap = binmap*newmap.fullMap
        newmap.fullMap =laplace(newmap.fullMap)
        return newmap

    def _surface_minimum_filter(self,c,window=3):
        # contour the map
        binmap = self.fullMap > float(c)
        # get the footprint array corresponding to 6 neighbors of a voxel
        fp = self._grid_footprint()
        # apply minimum filter to return surface voxels with zeros in
        binmap_surface = minimum_filter(binmap*1,footprint=fp,mode='constant',cval=0.0)
        # select those voxels with zeros filled in after applying minimum filter
        binmap_surface = ((binmap*1 - binmap_surface) == 1)*1

        return binmap_surface

    def _surface_features(self,c,window=21,itern=1):
        newmap = self.copy()
        binmap = self.fullMap > float(c)
        newmap.fullMap = binmap*1.0
        #print histogram(newmap.fullMap,30)
        for i in range(itern):
            newmap.fullMap = uniform_filter(newmap.fullMap,size=window,mode='constant',cval=0.0)
            newmap.fullMap = newmap.fullMap*binmap
            binmap = newmap.fullMap > 0.0
            minval = newmap.fullMap[binmap].min()
            newmap.fullMap = newmap.fullMap - minval + (0.001*minval)
            newmap.fullMap = newmap.fullMap*binmap
            #print len(newmap.fullMap[binmap])
            newmap.fullMap = newmap.fullMap/float(newmap.fullMap.max())
            #print histogram(newmap.fullMap,30)
        return newmap
    
    def _soft_mask(self,c=None,window=5,itern=3):
        #print histogram(newmap.fullMap,30)
        newmap = self.copy()
        if not c is None: newmap.fullMap = newmap.fullMap*(newmap.fullMap>float(c))
        #print '>>', c
        binmap = newmap.fullMap != 0 #only non-zero values replaced
        footprint_sph = self._make_spherical_footprint(window)
        #ma_map = ma.masked_greater_equal(newmap.fullMap,float(c)+newmap.fullMap.std())
        for i in range(itern):
            newmap.fullMap = generic_filter(newmap.fullMap,npmean,footprint=footprint_sph,mode='constant',cval=0.0)
            ##newmap.fullMap = uniform_filter(newmap.fullMap,size=window,mode='constant',cval=0.0)
        newmap.fullMap[binmap] = self.fullMap[binmap]
        return newmap.fullMap
    
    def _std_neigh(self,c=None,window=6):
        newmap = self.copy()
        if not c is None: newmap.fullMap = newmap.fullMap*(newmap.fullMap>float(c))
        footprint_sph = self._make_spherical_footprint(window)
        newmap.fullMap = generic_filter(newmap.fullMap,npstd,footprint=footprint_sph,mode='constant',cval=0.0)
        return newmap
    
    def _mean_neigh(self,c=None,window=6,itr=1):
        newmap = self.copy()
        if not c is None: newmap.fullMap = newmap.fullMap*(newmap.fullMap>float(c))
        footprint_sph = self._make_spherical_footprint(window)
        for i in range(itr):
            newmap.fullMap = generic_filter(newmap.fullMap,npmean,footprint=footprint_sph,mode='constant',cval=0.0)
        return newmap


    def _map_digitize(self,cutoff,nbins,left=False):
        try: from numpy import digitize
        except ImportError:
            print('Numpy Digitize missing, try v1.8')
        binMap = self.copy()
        bins = []
        step = (self.fullMap.max()-float(cutoff))/nbins
        ini = float(cutoff) + (0.000001*step)
        if left:
                ini = float(cutoff) - (0.000001*step)
        bins.append(ini)
        for ii in range(1,nbins+1):
                bins.append(float(cutoff) + ii*step)
        if bins[-1] < self.fullMap.max():
                bins = bins[:-1]
                bins.append(self.fullMap.max())

        for z in range(len(self.fullMap)):
            for y in range(len(self.fullMap[z])):
                binMap.fullMap[z][y] = digitize(self.fullMap[z][y],bins)
        return binMap

    def _map_depth(self,c):
        newmap = self.copy()
        binmap = self.fullMap > float(c)
        newmap.fullMap = binmap*1.0
        #print histogram(newmap.fullMap,30)
        for i in range(3):
                newmap.fullMap = uniform_filter(newmap.fullMap,size=21,mode='constant',cval=0.0)
                newmap.fullMap = newmap.fullMap*binmap
                binmap = newmap.fullMap > 0.0
                minval = newmap.fullMap[binmap].min()
                newmap.fullMap = newmap.fullMap - minval + 0.001
                newmap.fullMap = newmap.fullMap*binmap
                newmap.fullMap = newmap.fullMap/float(newmap.fullMap.max())
                #print histogram(newmap.fullMap,30)
        return newmap
    
        # remove smaller patches that can be found randomly w.r.t sizes from all patches
    def _label_patches(self,c,prob=0.2):
        fp = self._grid_footprint()
        binmap = self.fullMap > float(c)
        label_array, labels = measurements.label(self.fullMap*binmap,structure=fp)
        sizes = measurements.sum(binmap, label_array, range(labels + 1))
        if labels < 10:
            m_array = sizes < 0.05*sizes.max()
            ct_remove = numsum(m_array)
            remove_points = m_array[label_array]
            label_array[remove_points] = 0
            return (label_array>0)*(self.fullMap*binmap), labels-ct_remove+1

        means = measurements.mean(self.fullMap*binmap, label_array, range(1,labels + 1))
        freq,bins = histogram(sizes[1:],20)

        m_array = zeros(len(sizes))
        ct_remove = 0
        for i in range(len(freq)):
            fr = freq[i]
            s2 = bins[i+1]
            s1 = bins[i]
            p_size = float(fr)/float(numsum(freq))
            if s2 < 10 or p_size > prob:
                m_array = m_array+((sizes >= s1) & (sizes < s2))
                ct_remove += 1
        m_array = m_array>0
        remove_points = m_array[label_array]

        label_array[remove_points] = 0
        return (label_array>0)*(self.fullMap*binmap), labels-ct_remove
    
    # footprint array corresponding to 6 neighboring faces of a voxel
    def _grid_footprint(self):
        a = zeros((3,3,3))
        a[1,1,1] = 1
        a[0,1,1] = 1
        a[1,0,1] = 1
        a[1,1,0] = 1
        a[2,1,1] = 1
        a[1,2,1] = 1
        a[1,1,2] = 1

        return a
    #Added by AJ
    def _make_spherical_footprint(self,l):
        rad_z = arange(floor(l/2.0)*-1, ceil(l/2.0))
        rad_y = arange(floor(l/2.0)*-1, ceil(l/2.0))
        rad_x = arange(floor(l/2.0)*-1, ceil(l/2.0))

        rad_x = rad_x**2
        rad_y = rad_y**2
        rad_z = rad_z**2
        dist = srt(rad_z[:,None,None]+rad_y[:,None] + rad_x)
        #set_printoptions(threshold='nan')
        #print dist<=floor(l/2.0)
        return (dist<=floor(l/2.0))*1
    
    #Added by AJ
    # useful when a 'safe' contour level is chosen
    # for high resolution maps with grooves, this might erode part of useful density
    def _map_binary_opening(self,c,it=1):
        fp = self._grid_footprint()
        # current position can be updated based on neighbors only
        fp[1,1,1] = 0
        # threshold can be 1*std() to be safe?
        self.fullMap = self.fullMap*(self.fullMap > float(c))
        #print histogram(self.fullMap,10)
        self.fullMap = self.fullMap*binary_opening(self.fullMap>float(c),structure=fp,iterations=it)
        #print histogram(self.fullMap,10)


    def resize_map(self, new_size, centre=False):
        """
        
        Resize Map instance.
        
        
        Arguments:
           
           *new_size*
               3-tuple (x,y,z) giving the box size.
            *centre*
                default False
 
        Return:
            new Map instance
               
        """
        newMap = self.copy()
        newMap.fullMap = zeros(new_size)
        min_box = [min(x,y) for x,y in zip(newMap.box_size(), self.box_size())]
        newMap.fullMap[:min_box[0], :min_box[1], :min_box[2]] = self.fullMap[:min_box[0], :min_box[1], :min_box[2]]
        return newMap
    

    def normalise(self):
        """
        
        Return a new Map instance with normalised density values.
 
        Return:
            new Map instance         
       
        """
        newMap = self.copy()
        if self.fullMap.std() == 0:
            pass
        else:
            newMap.fullMap = (self.fullMap-self.fullMap.mean())/self.fullMap.std()
        return newMap
    
    #Added by AJ
    def pad_map(self,nx,ny,nz):
        """
        
        Pad a map (in place) with specified increments along each dimension.
        Arguments:
            *nx,ny,nz*
               Number of slices to pad in each dimension.
        Return:
            new Map instance         
       
        """
        gridshape = (self.fullMap.shape[0]+nz,self.fullMap.shape[1]+ny,self.fullMap.shape[2]+nx)
        new_array=zeros(gridshape)
        #min fill
        new_array.fill(self.fullMap.min())
        oldshape=self.fullMap.shape
        indz,indy,indx = round((gridshape[0]-oldshape[0])/2.),round((gridshape[1]-oldshape[1])/2.),round((gridshape[2]-oldshape[2])/2.)
        self.origin = (self.origin[0]-self.apix*indx,self.origin[1]-self.apix*indy,self.origin[2]-self.apix*indz)
        new_array[indz:indz+oldshape[0],indy:indy+oldshape[1],indx:indx+oldshape[2]] = self.fullMap
        self.fullMap = new_array

    def rotate_by_axis_angle(self, x, y, z, angle, CoM, rad=False):
        """
        
        Rotate the map around its centre given an axis and angle. 
        
        Arguments:
        
            *angle*
                angle (in radians if rad == True, else in degrees) to rotate map.
            *x,y,z*
                axis to rotate about, ie. x,y,z =  0,0,1 rotates the Map round the xy-plane.
            *CoM*
                centre of mass around which map will be rotated.
 
        Return:
            new Map instance         
               
        """
        
        m = Vector.axis_angle_to_matrix(x, y, z, angle, rad)
        # Calculate translation needed to rotate around CoM
        newCoM = CoM.matrix_transform(m.T)
        offset = CoM-newCoM
        # Apply transform
        newMap = self.matrix_transform(m, offset.x, offset.y, offset.z)
        return newMap

    def rotate_by_euler(self, x, y, z, CoM, rad=False):
        """
        
        Rotated map around pivot given by CoM using Euler angles. 
        
        Arguments:
        
            *x,y,z*
                Euler angles (in radians if rad == True, else in degrees) used to rotate map.
            *CoM*
                centre of mass around which map will be rotated.
            *x, y, z*
                translation in angstroms.
  
        Return:
            new Map instance         
          
        """
        m = Vector.euler_to_matrix(x, y, z, rad)
        # Calculate translation needed to rotate around CoM
        newCoM = CoM.matrix_transform(m.T)
        offset = CoM-newCoM
        # Apply transform
        newMap = self.matrix_transform(m, offset.x, offset.y, offset.z)
        return newMap
    
    #Added by AJ
    def _box_transform(self, mat):
        """
        Calculate box dimensions after rotation
        

        Arguments:

                *mat*
                        Input rotation matrix
        Return:
                new box shape
        """
        # Box corners
        v1 = Vector.Vector(self.origin[0],self.origin[1],self.origin[2])
        v2 = Vector.Vector(self.origin[0]+(self.apix*self.x_size()),self.origin[1],self.origin[2])
        v3 = Vector.Vector(self.origin[0]+(self.apix*self.x_size()),self.origin[1]+(self.apix*self.y_size()),self.origin[2])
        v4 = Vector.Vector(self.origin[0]+(self.apix*self.x_size()),self.origin[1]+(self.apix*self.y_size()),self.origin[2]+(self.apix*self.z_size()))
        v5 = Vector.Vector(self.origin[0],self.origin[1]+(self.apix*self.y_size()),self.origin[2])
        v6 = Vector.Vector(self.origin[0],self.origin[1],self.origin[2]+(self.apix*self.z_size()))
        v7 = Vector.Vector(self.origin[0],self.origin[1]+(self.apix*self.y_size()),self.origin[2]+(self.apix*self.z_size()))
        v8 = Vector.Vector(self.origin[0]+(self.apix*self.x_size()),self.origin[1],self.origin[2]+(self.apix*self.z_size()))

        # New corners
        v1 = v1.matrix_transform(mat)
        v2 = v2.matrix_transform(mat)
        v3 = v3.matrix_transform(mat)
        v4 = v4.matrix_transform(mat)
        v5 = v5.matrix_transform(mat)
        v6 = v6.matrix_transform(mat)
        v7 = v7.matrix_transform(mat)
        v8 = v8.matrix_transform(mat)

        output_shape=(self.x_size(),self.y_size(),self.z_size())
        max_x = 0
        max_y = 0
        max_z = 0
        ltmp = [v1,v2,v3,v4,v5,v6,v7,v8]
        # New ouput shape
        for i in range(8):
                for j in range(i,8):
                        if abs(ltmp[i].x - ltmp[j].x) > max_x:
                                max_x = abs(ltmp[i].x - ltmp[j].x)
                        if abs(ltmp[i].y - ltmp[j].y) > max_y:
                                max_y = abs(ltmp[i].y - ltmp[j].y)
                        if abs(ltmp[i].z - ltmp[j].z) > max_z:
                                max_z = abs(ltmp[i].z - ltmp[j].z)
        #output_shape = (int(max_x/self.apix),int(max_y/self.apix),int(max_z/self.apix))
        output_dimension = Vector.Vector(max_x,max_y,max_z)
        return output_dimension

    def _rotation_offset(self, mat, CoM1, CoM2, x=0,y=0,z=0,rad=False):
        newCoM = CoM2.matrix_transform(mat)
        offset = CoM1-newCoM
        return offset

    def rotate_by_matrix(self, mat, CoM, rad=False):
        """
        
        Rotated the map around pivot given by CoM using a rotation matrix 
        
        Arguments:
            *mat*
                3x3 matrix used to rotate map (in radians if rad == True, else in degrees).
            *CoM*
                rotation pivot, usually the centre of mass around which map will be rotated.
 
        Return:
            new Map instance         
       
        """
        # Calculate translation needed to rotate around CoM
        newCoM = CoM.matrix_transform(mat.T)
        offset = CoM-newCoM
        # Apply transform
        newMap = self.matrix_transform(mat, offset.x, offset.y, offset.z)
        return newMap
        
    def _matrix_transform_offset(self, mat, shape, x=0, y=0, z=0):
        newMap = self.copy()
        ## Transpose matrix
        #mat = mat.T
        # Convert offsets from angstrom to pixel values
        off_x = float(x/self.apix)
        off_y = float(y/self.apix)
        off_z = float(z/self.apix)
        newMap.fullMap = newMap.fullMap.swapaxes(0,2)
        newMap.fullMap = affine_transform(newMap.fullMap, mat, offset=(off_x, off_y, off_z), output_shape=shape, cval=(self.min()))
        newMap.fullMap = newMap.fullMap.swapaxes(0,2)
        return newMap


    def matrix_transform(self, mat, x=0, y=0, z=0):
        """
        
        Apply affine transform to the map.
                
        
        Arguments:
            
            *mat*
               affine 3x3 transformation matrix
            *shape*
            	new box dimensions
            *x, y, z*
                translation in angstroms.
            

        Return:
            new Map instance         
                
        """

        newMap = self.copy()
        # Transpose matrix
        mat = mat.T
        # Convert offsets from angstrom to pixel values
        off_x = x/self.apix
        off_y = y/self.apix
        off_z = z/self.apix
        # Calculate offsets due to rotation around (0,0,0) origin
        x_o = -self.x_origin()/self.apix
        y_o = -self.y_origin()/self.apix
        z_o = -self.z_origin()/self.apix

        off_x += x_o-mat[0,0]*x_o - mat[0,1]*y_o - mat[0,2]*z_o
        off_y += y_o-mat[1,0]*x_o - mat[1,1]*y_o - mat[1,2]*z_o
        off_z += z_o-mat[2,0]*x_o - mat[2,1]*y_o - mat[2,2]*z_o

        off_x = float(off_x)
        off_y = float(off_y)
        off_z = float(off_z)
        
        # Swap z and x axes, apply matrix, then swap back
        newMap.fullMap = newMap.fullMap.swapaxes(0,2)
        newMap.fullMap = affine_transform(newMap.fullMap, mat, offset=(off_x, off_y, off_z), cval=self.min())
        newMap.fullMap = newMap.fullMap.swapaxes(0,2)
        return newMap

    def change_origin(self, x_origin, y_origin, z_origin):
        """
        
        Change the origin of the map to a new origin. 
        
        Arguments:
            
        *x_origin, y_origin, z_origin*
            new co-ordinates of origin.

        Return:
            new Map instance         
        
        """
        newMap = self.copy()
        newMap.origin = (x_origin, y_origin, z_origin)
        return newMap

    def shift_origin(self, x_shift, y_shift, z_shift):
        """
        
        Shift the Map origin.

        Arguments:
            
        *x_origin, y_origin, z_origin*
            new co-ordinates of origin.
 
        Return:
            new Map instance         
       
        """

        newMap = self.copy()
        newMap.origin = (self.origin[0]+x_shift, self.origin[1]+y_shift, self.origin[2]+z_shift)
        return newMap

    def translate(self, x, y, z):
        """
        
        Translate the map by changing origin      
        
        Arguments:
            *x,y,z*
                 translation in angstroms
 
        Return:
            new Map instance         
       
        """
        sh = array([z/self.apix,y/self.apix,x/self.apix])
        newMap = self.copy()
        newMap.fullMap = shift(newMap.fullMap, sh, cval=self.min())
        #f_map = fftn(newMap.fullMap)
        #newMap.fullMap = real(ifftn((fourier_shift(f_map, shift))))
        return newMap
    
    def origin_change_maps(self,MapRef):
        """
        
        Return a new Map instance with origin changed accordingly to Reference Map 
            
        Arguments:
            *MapRef*
                Reference Map
        Return:
            new Map instance         
          
        """
        newMap = self.copy()
        origin_shift = [y-x for x,y in zip(newMap.origin, MapRef.origin)]
        #m2 = m2.shift_origin(origin_shift[0],origin_shift[1],origin_shift[2])
        newMap.translate(origin_shift[0],origin_shift[1],origin_shift[2])
        newMap.origin = MapRef.origin[:]
        return newMap
    

    def threshold_map(self, minDens, maxDens):
        """
        
        Create a Map instance containing only density values between the specified min and max values. 
        
        Arguments:
            
            *minDens*
                minimum density threshold
            *maxDens*
                maximum density threshold
  
        Return:
            new Map instance         
          
        """
        newMap1 = self.fullMap.copy()
        newMap1 = newMap1*(newMap1<maxDens)*(newMap1>minDens)
        newMap = self.copy()
        newMap.fullMap = newMap1
        return newMap
    
    #Added by AJ   
    def _find_level(self,vol):
        """
        Get the level corresponding to volume. 
        Arguments:
            *vol*
                volume within the contour
        Return:
            level corresponding to volume         
        """
        # initiate with reasonable values
        c1 = self.fullMap.min()
        vol_calc = float(vol)*2
        # loop until calc vol and exp vol matches
        it = 0
        flage = 0
        while ((vol_calc-float(vol))/(self.apix**3) > 10 and flage == 0):
            # mask only values >= previous sel
            mask_array = self.fullMap >= c1
            # compute histogram wt 1000 bins
            dens_freq,dens_bin = histogram(self.fullMap[mask_array],1000)
            # loop over bins to select a min level from the bins
            sum_freq = 0.0
            for i in range(len(dens_freq)):
                sum_freq += dens_freq[-(i+1)]
                dens_level = dens_bin[-(i+2)]
                vol_calc = sum_freq*(self.apix**3)
                # break when vol exceeds exp vol to select the bin level
                if vol_calc > float(vol):
                    sel_level = dens_level
                    bin_range = abs(dens_bin[-(i+1)]-dens_level)
                    it += 1
                    if sel_level <= c1: flage = 1
                    #print it, sel_level, c1, sum_freq, vol_calc, vol, flage
                    c1 = sel_level
                    if it == 3: flage = 1
                    break
        return sel_level
    
    #Added by AJ
    def _rotate_interpolate_to_grid(self,mat,gridshape,spacing_i,spacing_f,cort,fill=None,origin=None):
        """
        rotate a map and interpolate to new grid. 
        Arguments:
            *mat*
                rotation matrix
            *gridshape*
                new grid shape (z,y,x)
            *spacing_i*
                initial grid spacing
            *spacing_f*
                new grid spacing
            *cort*
                centre of rotation 
        Return:
            level corresponding to volume         
        """
        # grid to rotate and interpolate to : 
        # fill grid with min value
        # loop through and reverse rotate to find bounding coordinates
        nz = int(gridshape[0])
        ny = int(gridshape[1])
        nx = int(gridshape[2])
        map_centre = self.centre()
        if origin is None:
            origin = (map_centre.x - (nx*spacing_f)/2.0,map_centre.y - (ny*spacing_f)/2.0,map_centre.z - (nz*spacing_f)/2.0)
        orif = array(origin).view(float)
        orii = array(self.origin).view(float)
        nzi = int(self.fullMap.shape[0])
        nyi = int(self.fullMap.shape[1])
        nxi = int(self.fullMap.shape[2])
        si = float(spacing_i)
        sf = float(spacing_f)
        cor = array(cort).astype(float)
        new_grid = zeros(gridshape,dtype=float)
        if not fill == 0.0: new_grid.fill(self.min())
        maparray = self.fullMap.view(float)
        # rotation matrix transpose
        matt = (array(mat).T).astype(float)
        code = """
        /*calculate offset*/
        float corfx = ORIF1(0)+(nx*sf)/2.0;
        float corfy = ORIF1(1)+(ny*sf)/2.0;
        float corfz = ORIF1(2)+(nz*sf)/2.0;
        float crx = MATT2(0,0)*corfx + MATT2(0,1)*corfy + MATT2(0,2)*corfz;
        float cry = MATT2(1,0)*corfx + MATT2(1,1)*corfy + MATT2(1,2)*corfz;
        float crz = MATT2(2,0)*corfx + MATT2(2,1)*corfy + MATT2(2,2)*corfz;
        float offx = COR1(0) - crx;
        float offy = COR1(1) - cry;
        float offz = COR1(2) - crz;
        for (int z=0; z<nz; z++) {
          for (int y=0; y<ny; y++) {
            for (int x=0; x<nx; x++) {
              /*reverse rotation*/
              float xf = ORIF1(0)+(sf/2.0)+x*sf;
              float yf = ORIF1(1)+(sf/2.0)+y*sf;
              float zf = ORIF1(2)+(sf/2.0)+z*sf;
              float xi = MATT2(0,0)*xf + MATT2(0,1)*yf + MATT2(0,2)*zf;
              float yi = MATT2(1,0)*xf + MATT2(1,1)*yf + MATT2(1,2)*zf;
              float zi = MATT2(2,0)*xf + MATT2(2,1)*yf + MATT2(2,2)*zf;
              /*add offset*/
              xi = xi + offx;
              yi = yi + offy;
              zi = zi + offz;
              /*array coordinates*/
              xi = (xi - (ORII1(0)+(si/2.0)))/si;
              yi = (yi - (ORII1(1)+(si/2.0)))/si;
              zi = (zi - (ORII1(2)+(si/2.0)))/si;
              /*find bounds*/
              int x0 = floor(xi);
              int y0 = floor(yi);
              int z0 = floor(zi);
              int x1 = ceil(xi);
              int y1 = ceil(yi);
              int z1 = ceil(zi);
              /*std::cout << xf << ' ' << xi << ' ' << x0 << ' ' << x1 << ' ' << offx << std::endl;*/
              if ((x0 >= 0 && y0 >= 0 && z0 >= 0) && (x1 < nxi && y1 < nyi && z1 < nzi))
                { 
                float ma000 = MAPARRAY3(z0,y0,x0);
                float ma100 = MAPARRAY3(z1,y0,x0);
                float ma010 = MAPARRAY3(z0,y1,x0);
                float ma001 = MAPARRAY3(z0,y0,x1);
                float ma110 = MAPARRAY3(z1,y1,x0);
                float ma011 = MAPARRAY3(z0,y1,x1);
                float ma101 = MAPARRAY3(z1,y0,x1);
                float ma111 = MAPARRAY3(z1,y1,x1);
                NEW_GRID3(z,y,x) = ma000*(1-(xi-x0))*(1-(yi-y0))*(1-(zi-z0))+ma001*(xi-x0)*(1-(yi-y0))*(1-(zi-z0))+ma010*(1-(xi-x0))*(yi-y0)*(1-(zi-z0))+ma100*(1-(xi-x0))*(1-(yi-y0))*(zi-z0)+ma101*(xi-x0)*(1-(yi-y0))*(zi-z0)+ma110*(1-(xi-x0))*(yi-y0)*(zi-z0)+ma011*(xi-x0)*(yi-y0)*(1-(zi-z0))+ma111*(xi-x0)*(yi-y0)*(zi-z0);
                }
            }
          }
        }
        """
        weave.inline(code,['new_grid','matt','orii','orif','sf','si','nx','ny','nz','maparray','cor','nxi','nyi','nzi'], compiler = 'gcc',headers=["<math.h>"])
        self.fullMap = new_grid
        self.origin = origin
    
    #Added by AJ
    def _interpolate_to_grid(self,gridshape,s,ori,order_spline=3,fill='min'):
        """
        Spline interpolation to a grid.
        Arguments:
            
            *gridshape*
                shape of new grid array (z,y,x)
            *s*
                new grid spacing
            *ori*
                origin of the new grid
            *order_spine*
                order of the spline used for interpolation
        Return:
            Interpolated map
        """
        (ox,oy,oz) = (self.origin[0],self.origin[1],self.origin[2])
        (o1x,o1y,o1z) = (float(ori[0]),float(ori[1]),float(ori[2]))
        scale = float(s)/self.apix
        offset = (o1x-ox,o1y-oy,o1z-oz)

        new_map_origin = (o1x, o1y, o1z)
        grid_indices = indices(gridshape)
        z_ind = grid_indices[0]
        z_ind.ravel()
        y_ind = grid_indices[1]
        y_ind.ravel()
        x_ind = grid_indices[2]
        x_ind.ravel()
        z_ind = ((offset[2])/self.apix)+scale*z_ind
        y_ind = ((offset[1])/self.apix)+scale*y_ind
        x_ind = ((offset[0])/self.apix)+scale*x_ind
        #from datetime import datetime
        #print 'start', datetime.now().time()
        if order_spline > 1: filtered_array = spline_filter(self.fullMap,order=order_spline)
        else: filtered_array = self.fullMap
        if fill == 'zero': fillval = 0.0
        else: fillval = self.min()
        #filtered_array = self.fullMap
        new_array = map_coordinates(filtered_array,[z_ind,y_ind,x_ind],cval=fillval,order=order_spline,prefilter=False)
        new_map = Map(new_array.reshape(gridshape), new_map_origin, float(s), self.filename, self.header[:])

        #new_map.fullMap = new_array.reshape(gridshape)
        new_map.origin = (o1x, o1y, o1z)
        new_map.apix = float(s)
        return new_map
        '''
        # cubic smoothing polynome in 1D
        # Catmull-Rom style
        def interp_catmull_rom(p0, p1, p2, p3, f):
            return 0.5 * (p0 * (-f**3 + 2*f**2 - f) + p1 * (3*f**3 - 5*f**2 + 2) + p2 * (-3*f**3 + 4*f**2 + f) + p3 * (f**3 - f**2))
        '''
    
    #Added by AJ
    def _interpolate_to_grid1(self,gridshape,s,ori,order_spline=3,fill='min',bound=False):
        """
        Spline interpolation to a grid.
        Arguments:
            
            *gridshape*
                shape of new grid array (z,y,x)
            *s*
                new grid spacing
            *ori*
                origin of the new grid
            *order_spine*
                order of the spline used for interpolation
        Return:
            interpolated map
        """
        (ox,oy,oz) = (self.origin[0],self.origin[1],self.origin[2])
        (o1x,o1y,o1z) = (float(ori[0]),float(ori[1]),float(ori[2]))
        scale = float(s)/self.apix
        offset = (o1x-ox,o1y-oy,o1z-oz)

        new_map_origin = (o1x, o1y, o1z)
        '''
        grid_indices = indices(gridshape)
        z_ind = grid_indices[0]
        z_ind.ravel()
        y_ind = grid_indices[1]
        y_ind.ravel()
        x_ind = grid_indices[2]
        x_ind.ravel()
        '''
        # axis coordinates of new grid
        z_ind = arange(gridshape[0])
        y_ind = arange(gridshape[1])
        x_ind = arange(gridshape[2])
        # find coordinates relative to old grid
        z_ind = ((offset[2])/self.apix)+scale*z_ind
        y_ind = ((offset[1])/self.apix)+scale*y_ind
        x_ind = ((offset[0])/self.apix)+scale*x_ind
        # get indices of points inside the old grid
        z_mask_ind = ((nonzero((z_ind >= 0) & (z_ind < self.fullMap.shape[0]-1))[0])*1).astype(int)
        y_mask_ind = ((nonzero((y_ind >= 0) & (y_ind < self.fullMap.shape[1]-1))[0])*1).astype(int)
        x_mask_ind = ((nonzero((x_ind >= 0) & (x_ind < self.fullMap.shape[2]-1))[0])*1).astype(int)
        # indices of boundaries
        z_mask_ind1 = nonzero((z_ind >= self.fullMap.shape[0]-1)&(z_ind <  self.fullMap.shape[0]))[0]
        y_mask_ind1 = nonzero((y_ind >= self.fullMap.shape[1]-1)&(y_ind <  self.fullMap.shape[1]))[0]
        x_mask_ind1 = nonzero((x_ind >= self.fullMap.shape[2]-1)&(x_ind <  self.fullMap.shape[2]))[0]
        z_mask_ind0 = nonzero((z_ind < 0)&(z_ind > -1))[0]
        y_mask_ind0 = nonzero((y_ind < 0)&(y_ind > -1))[0]
        x_mask_ind0 = nonzero((x_ind < 0)&(x_ind > -1))[0]
        #searchsorted/floor
        # get the bound int coordinates 
        #searchsorted gives right bounds of orig array. substract 1 to get lower bound
        k_z = searchsorted(arange(self.fullMap.shape[0],dtype=float),z_ind-1,side='right')
        k_y = searchsorted(arange(self.fullMap.shape[1],dtype=float),y_ind-1,side='right')
        k_x = searchsorted(arange(self.fullMap.shape[2],dtype=float),x_ind-1,side='right')
        # new_grid 
        new_grid = zeros((gridshape))
        # extract lower bounds from original grid
        # check coordinate range
        x_grid1 = k_x[x_mask_ind].astype(int)
        #x_grid1[nonzero(x_grid1<0)] = 0
        y_grid1 = k_y[y_mask_ind].astype(int)
        #y_grid1[nonzero(y_grid1<0)] = 0
        z_grid1 = k_z[z_mask_ind].astype(int)
        #z_grid1[nonzero(z_grid1<0)] = 0
        # grid1 < size of orig map
        '''
        # meshgrid < 3D till numpy 1.6
        x_gridL,y_gridL,z_gridL = meshgrid(x_grid1,y_grid1,z_grid1)
        z_gridL.ravel()
        y_gridL.ravel()
        x_gridL.ravel()
        # OR use itertools
        from itertools import izip, product
        z_gridL,y_gridL,x_gridL = izip(*product(z_grid1,y_grid1,x_grid1))
        '''
        # faster option (array broadcasting - operations on different sizes)
        # indices from orig array
        tmp_grid = zeros((len(z_grid1),len(y_grid1),len(x_grid1)),dtype=int)
        z_gridL = (tmp_grid+z_grid1[...,newaxis,newaxis]).flatten()
        y_gridL = (tmp_grid+y_grid1[...,newaxis]).flatten()
        x_gridL = (tmp_grid+x_grid1).flatten()
        assert len(z_grid1) == len(z_mask_ind) and len(y_grid1) == len(y_mask_ind) and len(x_grid1) == len(x_mask_ind)
        # target array indices
        z_grid = (tmp_grid+z_mask_ind[...,newaxis,newaxis]).flatten()
        y_grid = (tmp_grid+y_mask_ind[...,newaxis]).flatten()
        x_grid = (tmp_grid+x_mask_ind).flatten()
        # fill with minimum value/zero
        if not fill == 'zero': new_grid.fill(self.min())
        # interpolate
        # have to check the following 3 lines : axis-wise
        '''
        # interpolate 1
        new_grid[z_grid,y_grid,x_grid] = (1.-(x_ind[x_grid]-x_gridL))*self.fullMap[z_gridL,y_gridL,x_gridL]+(x_ind[x_grid]-x_gridL)*self.fullMap[z_gridL,y_gridL,x_gridL+1]
        new_grid[z_grid,y_grid,x_grid] = (1-(y_ind[y_grid]-y_gridL))*new_grid[z_grid,y_grid,x_grid]+(y_ind[y_grid]-y_gridL)*((1.-(x_ind[x_grid]-x_gridL))*self.fullMap[z_gridL,y_gridL+1,x_gridL]+(x_ind[x_grid]-x_gridL)*self.fullMap[z_gridL,y_gridL+1,x_gridL+1])
        new_grid[z_grid,y_grid,x_grid] = (1-(z_ind[z_grid]-z_gridL))*new_grid[z_grid,y_grid,x_grid]
        new_grid[z_grid,y_grid,x_grid] = new_grid[z_grid,y_grid,x_grid] + (z_ind[z_grid]-z_gridL)*(((1.-(x_ind[x_grid]-x_gridL))*self.fullMap[z_gridL+1,y_gridL,x_gridL]+(x_ind[x_grid]-x_gridL)*self.fullMap[z_gridL+1,y_gridL,x_gridL+1])*(1-(y_ind[y_grid]-y_gridL))+((1.-(x_ind[x_grid]-x_gridL))*self.fullMap[z_gridL+1,y_gridL+1,x_gridL]+(x_ind[x_grid]-x_gridL)*self.fullMap[z_gridL+1,y_gridL+1,x_gridL+1])*(y_ind[y_grid]-y_gridL))
        '''
        from datetime import datetime
        # interpolation : C++
        nx = int(len(x_grid1))
        ny = int(len(y_grid1))
        nz = int(len(z_grid1))
        maparray = self.fullMap.view(float)
        xind = x_ind.view(float)
        yind = y_ind.view(float)
        zind = z_ind.view(float)
        code = """
        int k,j,i;
        int k1,j1,i1;
        for (int z=0; z<nz; z++) {
          for (int y=0; y<ny; y++) {
            for (int x=0; x<nx; x++) {
              k1 = Z_GRID11(z);
              j1 = Y_GRID11(y);
              i1 = X_GRID11(x);
              k = Z_MASK_IND1(z);
              j = Y_MASK_IND1(y); 
              i = X_MASK_IND1(x);
              float ma000 = MAPARRAY3(k1,j1,i1);
              float ma100 = MAPARRAY3(k1+1,j1,i1);
              float ma010 = MAPARRAY3(k1,j1+1,i1);
              float ma001 = MAPARRAY3(k1,j1,i1+1);
              float ma110 = MAPARRAY3(k1+1,j1+1,i1);
              float ma011 = MAPARRAY3(k1,j1+1,i1+1);
              float ma101 = MAPARRAY3(k1+1,j1,i1+1);
              float ma111 = MAPARRAY3(k1+1,j1+1,i1+1);
              float indx = XIND1(i);
              float indy = YIND1(j);
              float indz = ZIND1(k);
              NEW_GRID3(k,j,i) = ma000*(float) (1-(indx-i1))* (float) (1-(indy-j1))* (float) (1-(indz-k1)) +ma001*(indx-i1)*(1-(indy-j1))*(1-(indz-k1))+ma010*(1-(indx-i1))*(indy-j1)*(1-(indz-k1))+ma100*(1-(indx-i1))*(1-(indy-j1))*(indz-k1)+ma101*(indx-i1)*(1-(indy-j1))*(indz-k1)+ma110*(1-(indx-i1))*(indy-j1)*(indz-k1)+ma011*(indx-i1)*(indy-j1)*(1-(indz-k1))+ma111*(indx-i1)*(indy-j1)*(indz-k1);
            }
          }
        }"""
        # check 
        try:
            #print datetime.now().time() 
            weave.inline(code,['new_grid','x_grid1','y_grid1','z_grid1','z_mask_ind','y_mask_ind','x_mask_ind','xind','yind','zind','nx','ny','nz','maparray'],headers=["<math.h>"])
            #print datetime.now().time()
        except:
            print('C++ interpolation failed!, using alternate mode')
            # interpolation
            #print datetime.now().time()
            new_grid[z_grid,y_grid,x_grid] = (1.-(x_ind[x_grid]-x_gridL))*(1-(y_ind[y_grid]-y_gridL))*(1-(z_ind[z_grid]-z_gridL))*self.fullMap[z_gridL,y_gridL,x_gridL]+self.fullMap[z_gridL,y_gridL,x_gridL+1]*(x_ind[x_grid]-x_gridL)*(1-(y_ind[y_grid]-y_gridL))*(1-(z_ind[z_grid]-z_gridL))+self.fullMap[z_gridL,y_gridL+1,x_gridL]*(1.-(x_ind[x_grid]-x_gridL))*(y_ind[y_grid]-y_gridL)*(1-(z_ind[z_grid]-z_gridL))+self.fullMap[z_gridL+1,y_gridL,x_gridL]*(1.-(x_ind[x_grid]-x_gridL))*(1-(y_ind[y_grid]-y_gridL))*(z_ind[z_grid]-z_gridL)+self.fullMap[z_gridL+1,y_gridL,x_gridL+1]*(x_ind[x_grid]-x_gridL)*(1-(y_ind[y_grid]-y_gridL))*(z_ind[z_grid]-z_gridL)+self.fullMap[z_gridL+1,y_gridL+1,x_gridL]*(1.-(x_ind[x_grid]-x_gridL))*(y_ind[y_grid]-y_gridL)*(z_ind[z_grid]-z_gridL)+self.fullMap[z_gridL,y_gridL+1,x_gridL+1]*(x_ind[x_grid]-x_gridL)*(y_ind[y_grid]-y_gridL)*(1-(z_ind[z_grid]-z_gridL))+self.fullMap[z_gridL+1,y_gridL+1,x_gridL+1]*(x_ind[x_grid]-x_gridL)*(y_ind[y_grid]-y_gridL)*(z_ind[z_grid]-z_gridL)
            #print datetime.now().time()

        if bound == True:
            # note that the boundary interpolations are done just along one axis - to save time - not accurate
            # boundary interpolations
            for el in x_mask_ind1:
                new_grid[z_grid,y_grid,el] = (1-(x_ind[el]-k_x[el]))*self.fullMap[z_gridL,y_gridL,k_x[el]]+(x_ind[el]-k_x[el])*self.min()
            for el in y_mask_ind1:
                new_grid[z_grid,el,x_grid] = (1-(y_ind[el]-k_y[el]))*self.fullMap[z_gridL,k_y[el],x_gridL]+(y_ind[el]-k_y[el])*self.min()
            for el in z_mask_ind1:
                new_grid[el,y_grid,x_grid] = (1-(z_ind[el]-k_z[el]))*self.fullMap[k_z[el],y_gridL,x_gridL]+(z_ind[el]-k_z[el])*self.min()
            # 0's added on left extreme values outside the orig array
            k_x[x_mask_ind0] = -1.
            k_y[y_mask_ind0] = -1.
            k_z[z_mask_ind0] = -1.
            for el in x_mask_ind0:
                new_grid[z_grid,y_grid,el] = (1-(x_ind[el]-(-1)))*self.min()+(x_ind[el]-(-1))*self.fullMap[z_gridL,y_gridL,0]
            for el in y_mask_ind0:
                new_grid[z_grid,el,x_grid] = (1-(y_ind[el]-(-1)))*self.min()+(y_ind[el]-(-1))*self.fullMap[z_gridL,0,x_gridL]
            for el in z_mask_ind0:
                new_grid[el,y_grid,x_grid] = (1-(z_ind[el]-(-1)))*self.min()+(z_ind[el]-(-1))*self.fullMap[0,y_gridL,x_gridL]

        #interpolate
        new_map = Map(new_grid, new_map_origin, float(s), self.filename, self.header[:])

        #new_map.fullMap = new_array.reshape(gridshape)
        new_map.origin = (o1x, o1y, o1z)
        new_map.apix = float(s)
        return new_map
 
    #Added by AJ
    def _check_overlap(c1, c2, mat, cor, maxG):
        """
        Check whether transformation of a map overlaps with another.
        Note that for maps with large differences in x,y,z dimensions,
        some non-overlapping boxes may be returned but this is quick way
        to check for overlap without using the grids.
        
        Arguments:
            
            *c1*
                box centre of map1
            *c2*
                box centre of map2
            *mat*
                transformation matrix
            *cor*
                centre of rotation
            *maxG*
                maximum of the dimensions of the maps, in Angstrom
        Return:
            True/False
        """
        # centre of map2
        mat1 = matrix([[c2[0]],[c2[1]],[c2[2]]])
        # rotation matrix
        mat2 = matrix([mat[0][:-1],mat[1][:-1],mat[2][:-1]])
        c2_t = mat2*mat1
        # translate and adjust rotation offset
        c2_t[0] = c2_t[0]+mat[0][-1]+(float(cor[0])-c2_t[0])
        c2_t[1] = c2_t[1]+mat[1][-1]+(float(cor[1])-c2_t[1])
        c2_t[2] = c2_t[2]+mat[2][-1]+(float(cor[2])-c2_t[2])
        dist = math.sqrt(math.pow(c2_t[0]-c1[0],2)+math.pow(c2_t[1]-c1[1],2)+math.pow(c2_t[2]-c1[2],2))
        #print dist, maxG
        if dist < maxG/2.0: return True
        else: return False
    
    #Added by AJ
    def _mask_contour(self,thr,mar=0.0):
        '''
        Mask backgound beyond contour (in place)
        '''
        if mar != 0.0: level = thr-(mar*self.fullMap.std())
        else: level = thr
        minVal = self.min()
        a = ma.masked_less(self.fullMap, level, copy=True)
        self.fullMap = ma.filled(a, minVal)

    
    #Added by AJ,DV
    def _make_fourier_shell(self,fill=1.0):
        """
        For a given grid, make a grid with sampling frequencies as distance from centre
        Return:
          grid with sampling frequencies  
        """
        # numpy fftfreq : odd-> 0 to (n-1)/2 & -1 to -(n-1)/2 and even-> 0 to n/2-1 & -1 to -n/2 to check with eman
        # set frequencies in the range -0.5 to 0.5
        rad_z = arange(floor(self.z_size()/2.0)*-1, ceil(self.z_size()/2.0))/float(floor(self.z_size()))
        rad_y = arange(floor(self.y_size()/2.0)*-1, ceil(self.y_size()/2.0))/float(floor(self.y_size()))
        rad_x = arange(floor(self.x_size()/2.0)*-1, ceil(self.x_size()/2.0))/float(floor(self.x_size()))

        rad_x = rad_x**2
        rad_y = rad_y**2
        rad_z = rad_z**2

        dist = srt(rad_z[:,None,None]+rad_y[:,None] + rad_x)
        return dist

 
    #add by AP
    #for SS_CCF score
    # maybe change their name 
    def _get_maskArray(self, densValue):
        
        """ADDED by APP to use with SSCCC score"""
        mask_array = ma.masked_less_equal(self.fullMap, densValue)
        return ma.getmask(mask_array)
    
    
    def _get_maskMap(self, maskArray):
        
        """ADDED by APP to use with SSCCC score"""
        newMap= self.copy()
        newMap.fullMap *= 0
        masked_newMAP= ma.masked_array(self.fullMap, mask=maskArray,fill_value=0)
        #.filled() Return a copy of self, with masked values filled with a given value in this case 0
        newMap.fullMap=masked_newMAP.filled()
        return newMap
    
    def make_bin_map(self, cutoff):
        """
        
        Return a new Map instance that has been binarised. 
        All voxel with densities above and below the specified cutoff value are assigned a value of 1 and 0 respectively. 
        
        Arguments:
        
            *cutoff*
                cutoff density value

        Return:
            new binarised Map instance         
                
        """
        binMap = self.copy()
        #newMap = (self.fullMap > cutoff)
        binMap.fullMap = (self.fullMap > float(cutoff))*-1
        return binMap

    def _make_clash_map(self,apix=1.0):
        #note this is fine for the structure blured
        # look at note there DAVE NEED TO CTRL
##   that should be used in the blured but must be checked     
        """
        
        NOTE: NEEED TO BE CHECKED.
        
        Return an empty Map instance with set Angstrom per pixel sampling (default is 1)
 
        Return:
            new Map instance         
       
        """
        grid = zeros((self.box_size()[0]*self.apix/apix, self.box_size()[1]*self.apix/apix, self.box_size()[2]*self.apix/apix))
        clashMap = self.copy()
        clashMap.fullMap = grid
        clashMap.apix = apix
        return clashMap

    def resample_by_apix(self, new_apix):
        """
        
        Resample the map based on new_apix sampling.
        
        Arguments:
            *new_apix*
                Angstrom per pixel sampling
        
        Return:
            new Map instance         
     
        """
        new_map = self.copy()
        apix_ratio = self.apix/new_apix
        #new_map.apix = new_apix
        new_map.fullMap = resample(new_map.fullMap, self.z_size()*apix_ratio, axis=0)
        new_map.fullMap = resample(new_map.fullMap, self.y_size()*apix_ratio, axis=1)
        new_map.fullMap = resample(new_map.fullMap, self.x_size()*apix_ratio, axis=2)        
        new_map.apix = (self.apix*self.box_size()[2])/new_map.box_size()[2]
        return new_map

    def resample_by_box_size(self, new_box_size):
        """
        
        Resample the map based on new box size.
        
        Arguments
            *new_box_size*
                An array containing box dimension in ZYX format
                
        Return:
            new Map instance         
        
        """
        new_map = self.copy()
        new_map.fullMap = resample(new_map.fullMap, new_box_size[0], axis=0)
        new_map.fullMap = resample(new_map.fullMap, new_box_size[1], axis=1)
        new_map.fullMap = resample(new_map.fullMap, new_box_size[2], axis=2)        
        new_map.apix = (self.apix*self.box_size()[2])/new_box_size[2]
        return new_map

    # -- Map modifications involving filtering. All still return a new Map instance -- #

    def fourier_transform(self):
        """
        
        Apply Fourier transform on the density map.
        
        Return:
            new Map instance         
      
        """
        new_map = self.copy()
        new_map.fullMap = fftshift(fftn(self.fullMap))
        return new_map

    def laplace_filtered(self):
        """
        
       Apply Laplacian filter on density maps

        Return:
            new Map instance         
        
        """
        new_map = self.copy()
        new_map.fullMap = laplace(self.fullMap)
        return new_map
        

    # -- Methods that obtain information from the density map -- #

    def get_vectors(self):
        """
        
        Retrieve all non-zero density points in the form of Vector instances.
        
        Return:
            An array of 4-tuple (indices of the voxels in x,y,z format and density value)         
      
    """
        a = []
        for z in range(len(self.fullMap)):
            for y in range(len(self.fullMap[z])):
                for x in range(len(self.fullMap[z][y])):
                    if self.fullMap[z][y][x] != 0:
                        a.append((Vector.Vector((x*self.apix)+self.origin[0], (y*self.apix)+self.origin[1], (z*self.apix)+self.origin[2]),self.fullMap[z][y][x]))
        return array(a)
    
    def get_line_between_points(self, point1, point2):
        """
        
        Return an array of float values representing a line of density values between two points on the map.
        
        Arguments:
        
            *point1, point2* 
                Vector instances of the end points co-ordinates of the line.
 
        Return:
            array of floating values         
       """
        v1 = point1.minus(Vector.Vector(self.origin[0], self.origin[1], self.origin[2])).times(1.0/self.apix)
        v2 = point2.minus(Vector.Vector(self.origin[0], self.origin[1], self.origin[2])).times(1.0/self.apix)
        v3 = v2.minus(v1)
        noOfPoints = int(round(2*v3.mod()/self.apix))
        points = []
        for x in range(0, noOfPoints+1):
            p = v1.plus(v3.times(float(x)/noOfPoints))
            points.append(self.fullMap[p.z][p.y][p.x])
        return array(points)
    
    ### threshold map and find CoM
    def _get_com_threshold(self,c):
        """
            Return Vector instance of the centre of mass of the map.
            
        """

        newmap = self.copy()
        binmap = self.fullMap > float(c)
        newmap.fullMap = binmap*self.fullMap

        total_x_moment = 0.0
        total_y_moment = 0.0
        total_z_moment = 0.0
        total_mass = 0.0
        min_mass = newmap.min()
        #vectorMap = newmap.get_vectors()
        vectorMap = argwhere(newmap.fullMap)
        for v in vectorMap:
            #m = v[1]+min_mass
            m = newmap.fullMap[v[0]][v[1]][v[2]] + min_mass
            total_mass += m
            #total_x_moment += m*v[0].x
            #total_y_moment += m*v[0].y
            #total_z_moment += m*v[0].z
            total_x_moment += m*(self.origin[0]+v[2]*self.apix)
            total_y_moment += m*(self.origin[1]+v[1]*self.apix)
            total_z_moment += m*(self.origin[2]+v[0]*self.apix)
        x_CoM = total_x_moment/total_mass
        y_CoM = total_y_moment/total_mass
        z_CoM = total_z_moment/total_mass
        return Vector.Vector(x_CoM, y_CoM, z_CoM)


    def get_com(self):
        """
        
        Retrieve the centre of mass of the map.

        Return:
            Vector instance         
            
        """
        total_x_moment = 0.0
        total_y_moment = 0.0
        total_z_moment = 0.0
        total_mass = 0.0
        min_mass = self.min()
        vectorMap = self.get_vectors()
        for v in vectorMap:
            m = v[1]+min_mass
            total_mass += m
            total_x_moment += m*v[0].x
            total_y_moment += m*v[0].y
            total_z_moment += m*v[0].z
        x_CoM = total_x_moment/total_mass
        y_CoM = total_y_moment/total_mass
        z_CoM = total_z_moment/total_mass
        return Vector.Vector(x_CoM, y_CoM, z_CoM)

    def pixel_centre(self):
        """
        
        Return:
            Vector instance of the centre of the map in pixels.
       
        """
        x_centre = self.x_size()/2
        y_centre = self.y_size()/2
        z_centre = self.z_size()/2
        return Vector.Vector(x_centre, y_centre, z_centre)

    def centre(self):
        """
        Centre of the Map Instance
        
        Return:
            Vector instance of the centre of the map in Angstroms.
        
        """
        x_centre = self.origin[0]+(self.apix*self.x_size()/2)
        y_centre = self.origin[1]+(self.apix*self.y_size()/2)
        z_centre = self.origin[2]+(self.apix*self.z_size()/2)
        return Vector.Vector(x_centre, y_centre, z_centre)
    
    def mean(self):
        """
        
        Return:
            mean density value of map.
        
        """
        return self.fullMap.mean()

    def median(self):
        """
        
        Return:
            median density value of map.
        
        """
        return median(self.fullMap)

    def std(self):
        """
        
        Return:
            standard deviation of density values in map.
        
        """
        return self.fullMap.std()

    def min(self):
        """
        
        Return:
            minimum density value of the map.
        
        """
        return self.fullMap.min()
                    
    def max(self):
        """
        
        Return:
            maximum density value of the map.
        
        """
        return self.fullMap.max()

    def vectorise_point(self, x, y, z):
        """
        
        Return a tuple of the Angstrom co-ordinates and density value of a particular density point in map. 
        Transform the voxel specified by its indices (x,y,z) into a Vector object. The vector defines the position
        of the voxel with respect to the origin of the map. The magnitude of the vector is in Angstrom units.

        Arguments:
            *x, y, z*
                co-ordinates of the density point to be vectorised.
        
        Return:
            Vector instance
        """
        v_x = (self.apix*x)+self.origin[0]
        v_y = (self.apix*y)+self.origin[1]
        v_z = (self.apix*z)+self.origin[2]
        #density = self.fullMap[z][y][x]
        return Vector.Vector(v_x, v_y, v_z)

    def get_significant_points(self):
        """
        
        Retrieve all points with a density greater than one standard deviation above the mean.
        
        Return:
            An array of 4-tuple (indices of the voxels in x,y,z format and density value) 
        
        """
        sig_points = []
        boo = self.fullMap > (self.fullMap.mean() + self.fullMap.std())
        for z in range(self.z_size()):
            for y in range(self.y_size()):
                for x in range(self.x_size()):
                    if boo[z][y][x]:
                        sig_points.append(array([z,y,x,self.fullMap[z][y][x]]))
        return array(sig_points)

    def _get_random_significant_pairs(self, amount):
        """
        
        Return an array of tuple pairs of significant points randomly chosen from 'get_significant_points' function.
        For use with the DCCF and DLSF scoring functions.
        
        Arguments:
            *amount*
                number of significant point pairs to return.
        
        """
        sig_points = self.get_significant_points()
        sig_pairs = []
        size = len(sig_points)
        for r in range(amount):
            fst = sig_points[randrange(size)]
            snd = sig_points[randrange(size)]
            new_value = array([fst[0], fst[1], fst[2], snd[0], snd[1], snd[2], fst[3]-snd[3]])
            sig_pairs.append(new_value)
        return array(sig_pairs)

    def makeKDTree(self, minDens, maxDens):
        """
        
        Returns k-dimensional tree of points in the map with values between those chosen for quick nearest-neighbor lookup.
        
        Arguments:
            *minDens*
                minimum density value to include in k-dimensional tree.
            *maxDens*
               maximum density value to include in k-dimensional tree.
        
        Return:
            index into a set of k-dimensional points.
        """
        maplist = self.get_pos(minDens, maxDens)
        if len(maplist)!=0:
        	return KDTree(maplist)

    def get_pos(self, minDens, maxDens):
        """
        Identify a set of voxels in the map whose density values fall between the specified minimum and maximum values. 
                
        Arguments:
            *minDens*
                 minimum density value to include in array.
            *maxDens*
                maximum density value to include in array.
        
        Return:
            An array of 3-tuples (indices of the voxels in x,y,z format)
                
        """
        a = []
        for z in range(len(self.fullMap)):
            for y in range(len(self.fullMap[z])):
                for x in range(len(self.fullMap[z][y])):
                    if (self.fullMap[z][y][x] > minDens) and (self.fullMap[z][y][x] < maxDens):
                        a.append((x,y,z))
        return array(a)

    def _get_normal_vector(self, points):
        arr = points.view(int)
        vecnorm = zeros((len(points),3),dtype=float)
        np = len(points)
        xsize = int(self.x_size())
        ysize = int(self.y_size())
        zsize = int(self.z_size())
        mat = self.fullMap.view(float)
        
        #dims = nparray(ltmp,dtype=int)
        code = """
        int x,y,z;
        float mod;
        int flag=0;
        float xa,ya,za;
        for (int v=0; v<np; v++) {
            x = ARR2(v,2);
            y = ARR2(v,1);
            z = ARR2(v,0);
            /*std::cout << x << " " << y << " " << z << std::endl;*/
            xa = ya = za = 0.0;
            flag = 0;
            for (int i=x-1; i<x+2; i++) {
                for (int j=y-1; j<y+2; j++) {
                    for (int k=z-1; k<z+2; k++) {
                        if ((i == x) && (j == y) && (z == k))
                            continue;
                        else if ((i < 0) || (j < 0) || (k < 0))
                            continue;
                        else if ((i > xsize-1) || (j > ysize-1) || (k > zsize-1))
                            continue;
                        else {
                            if (MAT3(k,j,i) > MAT3(z,y,x)) {
                                /*std::cout << MAT3(z,y,x) << std::endl;*/
                                /*sum of unit vectors*/
                                mod = sqrt(pow(i-x,2)+pow(j-y,2)+pow(k-z,2));
                                if (mod != 0.0) {
                                    xa = xa + (i-x)/mod ;
                                    ya = ya + (j-y)/mod ;
                                    za = za + (k-z)/mod ;
                                    } 
                                flag = 1;
                                }
                            else if (MAT3(k,j,i) < MAT3(z,y,x))
                                flag = 1;
                            }
                    }
                }
            }
            if (flag != 0){
                /*unit average*/
                mod = sqrt(pow(xa,2)+pow(ya,2)+pow(za,2));
                if (mod != 0.0) {
                    VECNORM2(v,0) = xa/mod;
                    VECNORM2(v,1) = ya/mod;
                    VECNORM2(v,2) = za/mod;
                    }
                else {
                    VECNORM2(v,0) = 0.0;
                    VECNORM2(v,1) = 0.0;
                    VECNORM2(v,2) = 0.0;
                    }
                }
        }
           
        """
        # check
        try:
            #print datetime.now().time()
            weave.inline(code,['arr','np','mat','xsize','ysize','zsize','vecnorm'],headers=["<math.h>"])
            #print datetime.now().time()
            return vecnorm
        except:
            print('C++ NV scoring run failed!, using alternate mode')
            return None
        
    def get_normal_vector(self, x_pos, y_pos, z_pos):
        """
        
        Calculate the normal vector at the point specified. 
        Point calculated using 3SOM algorithm used by Ceulemans H. & Russell R.B. (2004).
       
        Arguments:
            *x_pos, y_pos, z_pos*
                pixel in map on which to calculate normal vector.
        Returns:
            Normal vector at the point specified
        """
	### check for no variationsi, if flag == 0 : no variation among neighbors
        flag = 0
        internal_vecs = []
        for x in range(x_pos-1, x_pos+2):
            for y in range(y_pos-1, y_pos+2):
                for z in range(z_pos-1, z_pos+2):
                    if (x_pos, y_pos, z_pos) == (x,y,z):
                        pass
                    elif x<0 or y<0 or z<0:
                        pass
                    elif x>self.x_size()-1 or y>self.y_size()-1 or z>self.z_size()-1:
                        pass
                    else:
                        if self.fullMap[z][y][x] > self.fullMap[z_pos][y_pos][x_pos]:
                            internal_vecs.append(Vector.Vector(x-x_pos,y-y_pos,z-z_pos).unit())
                            flag = 1
                        elif self.fullMap[z][y][x] < self.fullMap[z_pos][y_pos][x_pos]:
                            flag = 1

        sub_vector = Vector.Vector(0,0,0)
        for v in internal_vecs:
            sub_vector = sub_vector+v
        ### sometimes no vectors found (all zeroes )
        if len(internal_vecs) == 0 and flag == 0:
            return Vector.Vector(-9.0,-9.0,-9.0)
        return sub_vector.unit()

    def represent_normal_vectors(self, min_threshold, max_threshold):
        """
        
        Create a Structure instance representing normal vectors of density points specified.
                
        Arguments:
            *min_threshold, max_threshold*
                    minimum/maximum values to include in normal vector representation.
        Return:
            Structure Instance
        """
        atomList = []
        template = 'HETATM    1  C   NOR A   1      23.161  39.732 -25.038  1.00 10.00           C'
        m = self.copy()
        print(m.origin)
        for x in range(1, (m.box_size()[0]-1)):
            for y in range(1, (m.box_size()[1]-1)):
                for z in range(1, (m.box_size()[2]-1)):
                    if m.fullMap[z][y][x] > min_threshold and m.fullMap[z][y][x] < max_threshold:
                        #n_vec = m.get_normal_vector(x,y,z, min_threshold)
                        n_vec = m.get_normal_vector(x,y,z)
                        n_vec = n_vec.unit()
                        pos_vec = Vector.Vector((x*m.apix)+m.origin[0], (y*m.apix)+m.origin[1], (z*m.apix)+m.origin[2])
                        #a = BioPyAtom(template)
                        a = template.BioPyAtom()
                        a.x = pos_vec.x
                        a.y = pos_vec.y
                        a.z = pos_vec.z
                        b = BioPyAtom(template)
                        b.x = pos_vec.x + n_vec.x
                        b.y = pos_vec.y + n_vec.y
                        b.z = pos_vec.z + n_vec.z
                        c = BioPyAtom(template)
                        c.x = pos_vec.x + 2*n_vec.x
                        c.y = pos_vec.y + 2*n_vec.y
                        c.z = pos_vec.z + 2*n_vec.z
                        atomList.append(a)
                        atomList.append(b)
                        atomList.append(c)
        try:
            s = BioPy_Structure(atomList)
        except ZeroDivisionError:
            return atomList
        s.renumber_atoms()
        #for x in range(1, len(atomList),2):
        #    s.footer += 'CONECT   '+str(x)+'    '+str(x+1)+'\n'
        return s
    
    
#===============================================================================
#    
# #These two values should be calculated for the experimental map, and only
# need to be calculated once, at the beginning. Also note that when
# running the normal_vector_score, the experimental map should be the map1
# parameter (ie, the first one)     
#===============================================================================

    def get_point_map(self,min_thr,percentage=0.2):
        """
        
        Calculates the amount of point to use for the NV and CD score.
               
        Arguments:
  
            *min_thr*
                run get_primary_boundary on target map.
            *percentage*
                percentage of the protein volume.
        Return:
            
            number points
            
        """
        
        new_map = self.copy()
        prot_size = 1.*(new_map.fullMap > min_thr).sum()
        points = max(100, round(prot_size*percentage))
        return points


    @staticmethod
    def _gauss_bandpass(self,sigma,center=0.0):
        sigma = sigma/1.414 # for fourier space??
        map_filt = self.copy()
        size = self.x_size()
        if size <= self.y_size():
            size = self.y_size()
        if size <= self.z_size():
            size = self.z_size()
        flag_pyfftw = 1
        try:
            import pyfftw
        except ImportError: flag_pyfftw = 0
        #fft
        try:
            if flag_pyfftw == 0: raise ImportError
            inputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            outputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            fft = pyfftw.FFTW(inputa,outputa,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = self.fullMap[:,:,:]
            fft()
            ft = Map(fftshift(outputa), self.origin, self.apix, self.filename, self.header[:])
        except:
            ft = self.fourier_transform()
        #fourier shell, freq btw 0-0.5
        dist = self._make_fourier_shell(1.)
        # apply filter
        ft.fullMap = ft.fullMap*(np_exp(-((dist-center)**2)/(2*sigma*sigma)))
        #ifft
        try:
            if flag_pyfftw == 0: raise ImportError
            ifft = pyfftw.FFTW(inputa,outputa,direction='FFTW_BACKWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = ifftshift(ft.fullMap)[:,:,:]
            ifft()
            map_filt = Map(outputa.real.astype('float'), self.origin, self.apix, self.filename, self.header[:])
        except:
            map_filt.fullMap = real(ifftn(ifftshift(ft.fullMap)))
        return map_filt

    @staticmethod
    def _butterworth_lowpass(self,pass_freq,fall=0.0):
        eps = 0.882
        a = 10.624
        high_cutoff = 0.15*math.log10(1.0/pass_freq) + pass_freq # stop band frequency (used to determine the fall off)
        if fall == 0.0:
            fall = 2.0*(math.log10(eps/math.sqrt(a**2 - 1))/math.log10(pass_freq/float(high_cutoff)))
        cutoff_freq = float(pass_freq)/math.pow(eps,2/fall)
        map_filt = self.copy()
        size = self.x_size()
        if size <= self.y_size():
            size = self.y_size()
        if size <= self.z_size():
            size = self.z_size()
        flag_pyfftw = 1
        try:
            import pyfftw
        except ImportError: flag_pyfftw = 0
        #fft
        try:
            if flag_pyfftw == 0: raise ImportError
            inputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            outputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            fft = pyfftw.FFTW(inputa,outputa,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = self.fullMap[:,:,:]
            fft()
            ft = Map(fftshift(outputa), self.origin, self.apix, self.filename, self.header[:])
        except:
            ft = self.fourier_transform()
        # fourier shell, freq btw 0-0.5
        dist = self._make_fourier_shell(1.)
        # apply filter
        dist = dist/cutoff_freq
        dist = sqrt(1.0/(1.0+power(dist,fall)))
        ft.fullMap = ft.fullMap * dist
        #ifft
        try:
            if flag_pyfftw == 0: raise ImportError
            ifft = pyfftw.FFTW(inputa,outputa,direction='FFTW_BACKWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = ifftshift(ft.fullMap)[:,:,:]
            ifft()
            map_filt = Map(outputa.real.astype('float'), self.origin, self.apix, self.filename, self.header[:])
        except:
            map_filt.fullMap = real(ifftn(ifftshift(ft.fullMap)))
        return map_filt

    
    def _tanh_lowpass(self,cutoff,fall=0.2,ftmap=False):
        #cutoff = apix/reso is the stop band
        drop = math.pi/(2*float(cutoff)*float(fall))
        cutoff = min(float(cutoff),0.5)

        # fall determines smoothness of falloff, 0-> tophat, 1.0-> gaussian
        map_filt = self.copy()
        size = self.x_size()
        if size <= self.y_size():
            size = self.y_size()
        if size <= self.z_size():
            size = self.z_size()
        # if ft map is given
        if ftmap:
            #ft = self.copy()
            # fourier shells
            dist = self._make_fourier_shell(1.)
            # apply filter
            self.fullMap = self.fullMap * (0.5*(np_tanh(drop*(dist+cutoff))-np_tanh(drop*(dist-cutoff))))
            return

        flag_pyfftw = 1
        try:
            import pyfftw
        except ImportError:
            flag_pyfftw = 0

        try:
            if flag_pyfftw == 0: raise ImportError
            from datetime import datetime

            inputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            outputa = pyfftw.n_byte_align_empty(self.fullMap.shape, 16, 'complex128')
            fft = pyfftw.FFTW(inputa,outputa,direction='FFTW_FORWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = self.fullMap[:,:,:]
            #print 'pyfftw lowpass', datetime.now().time()
            fft()
            ft = Map(fftshift(outputa), self.origin, self.apix, self.filename, self.header[:])
          #print 'pyfftw lowpass end', datetime.now().time()
        except:
            ft = self.fourier_transform()
        # fourier shells
        dist = self._make_fourier_shell(1.)
        # apply filter
        ft.fullMap = ft.fullMap * (0.5*(np_tanh(drop*(dist+cutoff))-np_tanh(drop*(dist-cutoff))))
        '''
        # to view filter as a map
        tmpmap = self.copy()
        tmpmap.fullMap = (0.5*(np_tanh(drop*(dist+cutoff))-np_tanh(drop*(dist-cutoff))))
        ScoringFunctions_map.write_to_MRC_file(tmpmap,'filter_tanh.mrc')
        tmpmap.fullMap = dist
        ScoringFunctions_map.write_to_MRC_file(tmpmap,'radial.mrc')
        '''

        #ifft
        try:
            if flag_pyfftw == 0: raise ImportError
            ifft = pyfftw.FFTW(inputa,outputa,direction='FFTW_BACKWARD',axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
            inputa[:,:,:] = ifftshift(ft.fullMap)[:,:,:]
            ifft()
            map_filt = Map(outputa.real.astype('float'), self.origin, self.apix, self.filename, self.header[:])
        except:
            map_filt.fullMap = real(ifftn(ifftshift(ft.fullMap)))
        return map_filt

    @staticmethod
    def _tanh_bandpass(low_cutoff=0.0,high_cutoff=0.5,low_fall=0.1,high_fall=0.1):
        low_drop = math.pi/(2*float(high_cutoff-low_cutoff)*float(low_fall))
        high_drop = math.pi/(2*float(high_cutoff-low_cutoff)*float(high_fall))
        map_filt = self.copy()
        size = self.x_size()
        if size <= self.y_size():
            size = self.y_size()
        if size <= self.z_size():
            size = self.z_size()
        ft = self.fourier_transform()
        dist = self._make_fourier_shell(1.)
        ft.fullMap = ft.fullMap * 0.5*(np_tanh(high_drop*(dist+high_cutoff))-np_tanh(high_drop*(dist-high_cutoff))-np_tanh(low_drop*(dist+low_cutoff))+np_tanh(low_drop*(dist-low_cutoff)))
        map_filt.fullMap = real(ifftn(ifftshift(ft.fullMap)))
        return map_filt


    def get_primary_boundary(self, molWeight, low, high,vol_factor=1.21):
        """
        
        Calculates primary boundary density value. Volume of pixels with greater density than output is 
        equivalent to volume given by molecular weight of protein. Uses recursive algorithm. 
        NOTE: when used to calculated the NV score this value should be calculated for the experimental map.
               
        Arguments:
  
            *molWeight*
                molecular weight of protein;
                use get_prot_mass_from_atoms() if your structure contains HETATOMS else use get_prot_mass_from_res().
           *low, high*
               minimum and maximum values between which the boundary will be taken. 
               Initial values should be given by minimum and maximum density values in map.
            *vol_factor*
                in cubic Angstroms per Dalton. 
                This is the approximate value for globular proteins used in Chimera (Petterson et al, 2004) from Harpaz 1994.
                Other recommended volume factor are 1.5 (1.1-1.9) cubic Angstroms per Dalton in EMAN Volume/mass conversions assume a density of 1.35 g/ml (0.81 Da/A3) (~1.23A3/Da)
        
        Return:
            
            primary boundary density value (float)
            
        """
        #print low,high

        new_map = self.copy()
        if numsum(new_map.fullMap == low) > new_map.fullMap.size/10:
            all_dens = new_map.fullMap.flatten()
            all_dens = set(all_dens)
            all_dens=sorted(all_dens)
            l_ind = all_dens.index(low)
            low = all_dens[l_ind+1]
        if high-low < 0.0000002 or high < low:
            est_weight = long(numsum(new_map.fullMap > low)*new_map.apix**3/(vol_factor*1000))
            print('Exact molecular weight cannot be found. Approx. weight of '+str(est_weight)+' used instead.') 
            
            return low
        thr = low+(high-low)/2
        this_weight = long(numsum(new_map.fullMap > thr)*new_map.apix**3/(vol_factor*1000))
        #print "this_weight",this_weight, thr
        #this_weight = long((self.fullMap > thr).sum()*self.apix**3/1210)
        #print thr, this_weight, this_weight.sum()
        if this_weight == long(molWeight):
            return thr
        elif this_weight > molWeight:
            return new_map.get_primary_boundary(molWeight, thr, high)
        elif this_weight < molWeight:
            return new_map.get_primary_boundary(molWeight, low, thr)

    def _get_second_boundary_outward(self, primary_boundary, noOfPoints, low, high, err_percent=1):
        
        """
        PRIVATE FUNCTION to calculate the second bound density value. 
        Searching from primary boundary outward.
        For a given boundary value, it calculates the second bound density value 
        such that a specified number of points whose density values fall between the defined boundaries
        Uses recursive algorithm. 
        Arguments:  
           *primary_boundary*
               primary threshold, normally given by get_primary_boundary method based on protein molecular weight.
           *noOfPoints*
                   Number of points to use in the normal vector score - try first with 10% (maybe 5%better) of the number of points in the map ( round((self.map_size())*0.1)
           *low, high*
               minimum and maximum values between which the threshold will be taken.
               low should be equal to the value returned by the get_primary_boundary() method and high is the maximum density values in map.  
           *err_percent*
                 default value of 1. Allowed to find a secondary boundary that includes a 1% error.                
        
        Return:
             outward secondary boundary density value           
        
        """
          
        #print low, high
        if high-low < 0.0000002 or high<low:
            est_weight =  numsum((self.fullMap < low)*(self.fullMap > primary_boundary))
            print('Non optimal number of pixels to match. Try changing this value or increasing the size of err_percent ')
            return 1j
        
        thr = low+(high-low)/2
        this_no_of_points = numsum((self.fullMap < thr)*(self.fullMap > primary_boundary))
        #print thr, this_no_of_points
        if  abs(this_no_of_points - noOfPoints) < err_percent*noOfPoints/100.:
            return thr
        elif this_no_of_points < noOfPoints:
            return self._get_second_boundary_outward(primary_boundary, noOfPoints, thr, high)
        elif this_no_of_points > noOfPoints:
            return self._get_second_boundary_outward(primary_boundary, noOfPoints, low, thr)

    def _get_second_boundary_inward(self, primary_boundary, noOfPoints, low, high, err_percent=1):
        
        """
        
        PRIVATE FUNCTION to calculate the second bound density value. 
        Searching from primary boundary inward.
        For a given boundary value, it calculates the second bound density value 
        such that a specified number of points whose density values fall between the defined boundaries
        Uses recursive algorithm. 
        Arguments:  
           *primary_boundary*
               primary threshold, normally given by get_primary_boundary method based on protein molecular weight.
           *noOfPoints*
                   Number of points to use in the normal vector score - try first with 10% (maybe 5%better) of the number of points in the map ( round((self.map_size())*0.1)
           *low, high*
               minimum and maximum values between which the threshold will be taken.
               low should be equal to the value returned by the get_primary_boundary() method and high is the maximum density values in map.  
           *err_percent*
                 default value of 1. Allowed to find a secondary boundary that includes a 1% error.
        
        Return:
             outward secondary boundary density value                       
            """
                       
#         print low, high
        if high-low < 0.0000002 or high<low:
            est_weight =  numsum((self.fullMap < low)*(self.fullMap > primary_boundary))
            print('Non optimal number of pixels to match. Try changing this value or increasing the size of err_percent ')
            return 1j
        
        thr = high-(high-low)/2
        this_no_of_points = numsum((self.fullMap < thr)*(self.fullMap > primary_boundary))
        #print thr, this_no_of_points, noOfPoints
        if  abs(this_no_of_points - noOfPoints) < err_percent*noOfPoints/100.:
            return thr
        elif this_no_of_points < noOfPoints:
            return self._get_second_boundary_inward(primary_boundary, noOfPoints, thr, high)
        elif this_no_of_points > noOfPoints:
            #print 'test'
            return self._get_second_boundary_inward(primary_boundary, noOfPoints, low, thr)

    def get_second_boundary(self, primary_boundary, noOfPoints, low, high, err_percent=1):
        
        """
        Calculate the second bound density value. For a given boundary value, 
        it calculates the second bound density value such that a specified number
        of points whose density values fall between the defined boundaries
        Uses recursive algorithm. 
         
        Arguments:  
           *primary_boundary*
               primary threshold, normally given by get_primary_boundary method based on protein molecular weight.
           *noOfPoints*
                   Number of points to use in the normal vector score - try first with 10% (maybe 5%better) of the number of points in the map ( round((self.map_size())*0.1)
           *low, high*
               minimum and maximum values between which the threshold will be taken.
               low should be equal to the value returned by the get_primary_boundary() method and high is the maximum density values in map.  
           *err_percent*
                 default value of 1. Allowed to find a secondary boundary that includes a 1% error.
        
        Return:
            secondary boundary density value           
        """
        
        bou = self._get_second_boundary_outward(primary_boundary, noOfPoints, low, high, err_percent)
        if bou == 1j:
            bou = self._get_second_boundary_inward(primary_boundary, noOfPoints, low, high, err_percent)
        return bou
    
    # DO THIS (SHRINK BOX SIZE TO AROUND SD VALUE)
    def _shrink_map(self, sd=2.):
        pass

    # -- Writing Map instance to files -- #
    
    def _write_to_xplor_file(self, xplorFileName):
        """OBSOLETE PRIVATE FUNCTION
        
        xplorFileName = name of file to write to. Note that this function does not automatically append a .xplor suffix."""
        xplor = '\n 1 !NTITLE\n'
        
        xplor += 'REMARKS '+'"' + xplorFileName + '"' + '    written by ME!\n'
        
        xplor += self._pad_grid_line_no(self.z_size()) + self._pad_grid_line_no(0) + self._pad_grid_line_no(self.z_size()-1) + \
                 self._pad_grid_line_no(self.y_size()) + self._pad_grid_line_no(0) + self._pad_grid_line_no(self.y_size()-1) + \
                 self._pad_grid_line_no(self.x_size()) + self._pad_grid_line_no(0) + self._pad_grid_line_no(self.x_size()-1) + '\n'
        
        xplor += self._convert_point_to_string(self.apix*self.z_size()) + self._convert_point_to_string(self.apix*self.y_size()) + \
                 self._convert_point_to_string(self.apix*self.x_size())
        
        xplor += self._convert_point_to_string(90.0) + self._convert_point_to_string(90.0) + self._convert_point_to_string(90.0) + '\n'
        xplor += 'ZYX\n'
        flatList = self.fullMap.flatten()
        blockSize = self.z_size()*self.y_size()
        blockNo = 0
        offset = 0
        for point in range(len(flatList)):
            if ((point-offset)%6 == 0 and point%blockSize != 0): 
                    xplor += '\n'
            if point%blockSize == 0:
                xplor += '\n'+ self._pad_grid_line_no(blockNo) +'\n'
                blockNo += 1
                offset = point%6
            xplor += self._convert_point_to_string(real(flatList[point]))
        xplor += '\n   -9999'
        f = file(xplorFileName, 'w')
        f.write(xplor)
        f.close()

    def _write_to_situs_file(self, situsFileName):
        """One day I will do this."""
        pass

    def write_to_MRC_file(self, mrcfilename):
        
        """
        
        Write out a MRC file
        
        Arguments:
            *mrcfilename*
                name of the output mrc file

        """
        h = self.header_to_binary()
        maparray = array(self.fullMap, dtype='float32')
        f = open(mrcfilename, 'wb')
        f.write(h)
        f.write(maparray.tostring())
        f.close()

    def update_header(self):
        """
        
        Update self.header to values currently relevant.
        
        """
        nx = int32(self.x_size())
        ny = int32(self.y_size())
        nz = int32(self.z_size())
        mode = int32(2)


        #nxstart = int32(self.origin[0]/self.apix)
        #nystart = int32(self.origin[1]/self.apix)
        #nzstart = int32(self.origin[2]/self.apix)
        
        
        nxstart = int32(-nx/2)
        nystart = int32(-ny/2)
        nzstart = int32(-nz/2)
        mx = nx
        my = ny
        mz = nz
        xlen = float32(nx*self.apix)
        ylen = float32(ny*self.apix)
        zlen = float32(nz*self.apix)
        alpha = float32(90)
        beta = float32(90)
        gamma = float32(90)
        mapc = int32(1)
        mapr = int32(2)
        maps = int32(3)
        amin = float32(self.min())
        amax = float32(self.max())
        amean = float32(self.mean())
        ispg = int32(0)
        nsymbt = int32(0)
        extra = ''
        xorigin = float32(self.origin[0])
        yorigin = float32(self.origin[1])
        zorigin = float32(self.origin[2])
        mapword = 'MAP '
        if sys.byteorder == 'little':
            byteorder = 0x44440000
        else:
            byteorder = 0x11110000
        rms = float32(self.std())
        nlabels = int32(1)
        label0 = 'Created by TEMpy on: ' + str(datetime.date.today())
        otherlabels = ''
        
        self.header = [nx, ny, nz, mode, nxstart, nystart, nzstart, mx, my, mz, xlen, ylen, zlen, alpha, beta, gamma,\
                       mapc, mapr, maps, amin, amax, amean, ispg, nsymbt, extra, xorigin, yorigin, zorigin, mapword, byteorder,\
                       rms, nlabels, label0, otherlabels]


    def header_to_binary(self):
        """
        
        Returns binary version of map header data. For use in writing out density maps in MRC file format. 
        
        """
        nx = int32(self.x_size())
        ny = int32(self.y_size())
        nz = int32(self.z_size())
        mode = int32(2)
        #nxstart = int32(-nx/2)
        #nystart = int32(-ny/2)
        #nzstart = int32(-nz/2)
        nxstart = int32(self.origin[0]/self.apix)
        nystart = int32(self.origin[1]/self.apix)
        nzstart = int32(self.origin[2]/self.apix)
        
        mx = nx
        my = ny
        mz = nz
        xlen = float32(nx*self.apix)
        ylen = float32(ny*self.apix)
        zlen = float32(nz*self.apix)
        alpha = int32(90)
        beta = int32(90)
        gamma = int32(90)
        mapc = int32(1)
        mapr = int32(2)
        maps = int32(3)
        amin = float32(self.min())
        amax = float32(self.max())
        amean = float32(self.mean())
        ispg = int32(0)
        nsymbt = int32(0)
        extra = zeros(10).tostring()
        xorigin = float32(self.origin[0])
        yorigin = float32(self.origin[1])
        zorigin = float32(self.origin[2])
        mapword = 'MAP '
        if sys.byteorder == 'little':
            byteorder = 0x44440000
        else:
            byteorder = 0x11110000
        rms = float32(self.std())
        nlabels = int32(1)
        label0 = 'Created by TEMpy on: ' + str(datetime.date.today()) + zeros(49).tostring()
        otherlabels = zeros(90).tostring()

        fm_string = '=10l6f3l3f2l100s3f4slfl80s720s'
        packed = binary.pack(fm_string, nx, ny, nz, mode, nxstart, nystart, nzstart, mx, my, mz, xlen, ylen, zlen, alpha, beta, gamma,\
                       mapc, mapr, maps, amin, amax, amean, ispg, nsymbt, extra, xorigin, yorigin, zorigin, mapword, byteorder,\
                       rms, nlabels, label0, otherlabels)
        return packed

    def _pad_grid_line_no(self, no):
        """Private function to help write data to map files."""
        s = str(no)
        spaces = ''
        for x in range(8-len(s)):
            spaces += ' '
        s = spaces + s
        return s
    
    def _convert_point_to_string(self, point):
        """Private function to help write data to map files."""
        exp = 0
        sign = ''
        if(abs(point)<0.0001):
            point = 0.0
        if point >= 0:
            sign = '+'
        else:
            sign = '-'
        while abs(point) >= 10.0:
            point /= 10.0
            exp += 1
        pointString = str(point)
        if len(pointString) < 7:
            for x in range(len(pointString), 7):
                pointString += '0'
        elif len(pointString) > 7:
            pointString = pointString[:7]
        pointString += 'E' + sign + '0' + str(exp)
        return ' '+pointString

    def _get_component_volumes(self,struct, apix, blurrer):
        """
        Private function for check on Map instance.
        Return a a list containing the volume of the individual components based on a grid with voxel size set to apix
                        
        Arguments:  
            *struct*
                Structure object containing one or more chains.
            *apix*
                voxel size of the grid.
            *blurrer*
                Instance of a StructureBlurrer class. 
                
        Return:
            Return a a list containing the volume of the individual components         

        """
        mapCoM = self.get_com()
        ssplit = struct.split_into_chains()
        temp_grid = self._make_clash_map(apix)

        overlay_maplist = []
        cvol = []
        for x in ssplit:
            tx = mapCoM[0] - x.CoM[0]
            ty = mapCoM[1] - x.CoM[1]
            tz = mapCoM[2] - x.CoM[2]
            x.translate(tx,ty,tz)
            overlay_maplist.append(blurrer.make_atom_overlay_map1(temp_grid, x))
        for y in overlay_maplist:
        	cvol.append(y.fullMap.sum()*(apix**3))
        return cvol
 
 #added by PAP
    def map_rotate_by_axis_angle(self, x, y, z, angle, CoM, rad=False):
        """
        Return a new Map instance rotated around its centre.
                        
        Arguments:  
            *angle*
                angle (in radians if rad == True, else in degrees) to rotate map.
            *x,y,z*
                axis to rotate about, ie. x,y,z =  0,0,1 rotates the Map round the xy-plane.
            *CoM*
                centre of mass around which map will be rotated. 
                
        Return:
            Rotated new Map instance         
        """

        m = Vector.axis_angle_to_matrix(x, y, z, angle, rad)
        #print Vector(x,y,z).unit()

        # Calculate translation needed to rotate around CoM
        newCoM = CoM.matrix_transform(m.T)
        offset = CoM-newCoM

        # Apply transform
        newMap = self.matrix_transform(m, offset.x, offset.y, offset.z)
        return newMap
