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


from numpy import array,  zeros, real,sqrt,exp, mgrid, transpose
from scipy.fftpack import fftn, ifftn
from scipy.ndimage import fourier_gaussian,gaussian_filter,uniform_filter
from TEMPy.EMMap import Map

class StructureBlurrer:
    """ 
    
    A class to generates a density map from a structure instance.
    
    """

    def __init__(self):
        pass

    def protMap(self, struct, apix, resolution=None,filename="None"):
        
        """
        
        Returns an Map instance sized and centred based on the atomic structure.
        
        Arguments:
        
           *apix*
               Angstroms per pixel for the Map to be outputted.
           *resolution*
                Target resolution of the outputted map.
           *sigma_coeff*
               Sigma width of the Gaussian used to blur the atomic structure.
           *filename* 
               output name of the map file.
               
           """

        # Build empty template map based on the size of the protein and the resolution.
        extr = struct.get_extreme_values()
        if not resolution is None: edge = int(2*resolution/apix)+4
        else: edge = 10
        x_size = int((extr[1]-extr[0])/apix)+edge
        y_size = int((extr[3]-extr[2])/apix)+edge
        z_size = int((extr[5]-extr[4])/apix)+edge

        # Origin calculated such that the centre of the map is the centre of mass of the protein.
        half_x = max(struct.CoM.x - extr[0],extr[1]-struct.CoM.x)
        ##x_origin = struct.CoM.x-(apix*x_size/2.0)
        if half_x < (apix*x_size/2.0): half_x = apix*x_size/2.0
        x_origin = struct.CoM.x - half_x - edge*apix
        # apj: if com is not near the geometric centre of protein
        x_size = int(half_x*2.0/apix + 2*edge)
        ##y_origin = struct.CoM.y-(apix*y_size/2.0)
        # apj: if com is not near the geometric centre of protein
        half_y = max(struct.CoM.y - extr[2],extr[3]-struct.CoM.y)
        if half_y < (apix*y_size/2.0): half_y = (apix*y_size/2.0)
        y_origin = struct.CoM.y - half_y - edge*apix
        y_size = int(half_y*2.0/apix+ 2*edge)
        ##z_origin = struct.CoM.z-(apix*z_size/2.0)
        # apj: if com is not near the geometric centre of protein
        half_z = max(struct.CoM.z - extr[4],extr[5]-struct.CoM.z)
        if half_z < (apix*z_size/2.0): half_z = apix*z_size/2.0
        z_origin = struct.CoM.z - half_z - edge*apix
        z_size = int(half_z*2.0/apix+ 2*edge)
        
        newMap = zeros((z_size, y_size, x_size))
        fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
        return fullMap
      #add by IF
    
    def protMapBox(self, struct, apix, resolution,box_size_x,box_size_y,box_size_z,filename):
        """
        Create a Map instance sized and centered based on the atomic structure.
        
        
        Arguments:
        
            *struct*
                the Structure instance.
            *apix*
                Angstroms per pixel for the output Map.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                x dimension of output map box in Angstroms.
            *box_size_y*
                y dimension of output map box in Angstroms.
            *box_size_z*
                z dimension of output map box in Angstroms.
            *filename*
                output name of the map file.
        
        Return:
            A Map instance
            
        """

        # Build empty template map based on the size of the protein and the resolution.
        x_size = int(box_size_x)
        y_size = int(box_size_y)
        z_size = int(box_size_z)

        # Origin calculated such that the centre of the map is the centre of mass of the protein.
        x_origin = struct.CoM.x-(apix*x_size/2.0)
        y_origin = struct.CoM.y-(apix*y_size/2.0)
        z_origin = struct.CoM.z-(apix*z_size/2.0)
        
        newMap = zeros((z_size, y_size, x_size))
        fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
        return fullMap
  
    def mapGridPosition(self, densMap, atom):
        
        """
        
        Returns the index of the nearest pixel to an atom, and atom mass (4 values in list form).
        
        Arguments:
        
           *densMap*
               Map instance the atom is to be placed on.
           *atom*
               Atom instance.
               
           """
        origin = densMap.origin
        apix = densMap.apix
        box_size = densMap.box_size()
        x_pos = int(round((atom.x-origin[0])/apix,0))
        y_pos = int(round((atom.y-origin[1])/apix,0))
        z_pos = int(round((atom.z-origin[2])/apix,0))
        #print "grid_pos", x_pos,y_pos,z_pos,atom.x-origin[0], atom.y-origin[1], atom.z-origin[2]
        
        #MODIFIED BY PAP
        if((box_size[2] > x_pos >= 0) and (box_size[1] > y_pos >= 0) and (box_size[0] > z_pos >= 0)):
            return (x_pos, y_pos, z_pos, atom.mass)
        else:
            return 0
    
    def mapGridPositions_vdw(self, densMap, atom, gridtree):
        
        """
        
        Returns the index of the nearest pixel to an atom, and atom mass (4 values in list form).
        
        Arguments:
        
           *densMap*
               Map instance the atom is to be placed on.
           *atom*
               Atom instance.
               
           """
        origin = densMap.origin
        apix = densMap.apix
        box_size = densMap.box_size()
        x_pos = int(round((atom.x-origin[0])/apix,0))
        y_pos = int(round((atom.y-origin[1])/apix,0))
        z_pos = int(round((atom.z-origin[2])/apix,0))
        if((densMap.x_size() > x_pos >= 0) and (densMap.y_size() > y_pos >= 0) and (densMap.z_size() > z_pos >= 0)):
            #search all points withing 1.5sigma
            list_points = gridtree.query_ball_point([atom.x,atom.y,atom.z], atom.vdw)
            return list_points, (x_pos, y_pos, z_pos)
        else:
            print("Warning, atom out of map box")
            return [], ()

    def maptree(self,densMap,strmap=None):
        """
        
        Returns the KDTree of coordinates from a map grid.
        
        Arguments:
        
           *densMap*
               Map instance the atom is to be placed on.
               
        """
        origin = densMap.origin
        apix = densMap.apix
        box_size = densMap.box_size()
        nz,ny,nx = densMap.fullMap.shape

        #convert to real coordinates 
        zg,yg,xg = mgrid[0:nz,0:ny,0:nx]
        #to get indices in real coordinates
        zg = zg*apix + origin[2] + apix/2.0
        yg = yg*apix + origin[1] + apix/2.0
        xg = xg*apix + origin[0] + apix/2.0
        indi = zip(xg.ravel(), yg.ravel(), zg.ravel())
        try: 
            from scipy.spatial import cKDTree
            gridtree = cKDTree(indi)
        except ImportError:
            try:
                from scipy.spatial import KDTree
                gridtree = KDTree(indi)
            except ImportError: return
        return gridtree


    #voxels occupied by the atom density (3 sigma gaussian)
    def mapGridPositions(self,densMap, atom, gridtree,res_map,sim_sigma_coeff=0.187):
        """
        
        Returns the indices of the nearest pixels to an atom as a list.
        
        Arguments:
        
           *densMap*
               Map instance the atom is to be placed on.
           *atom*
               Atom instance.
           *gridtree*
               KDTree of the map coordinates (absolute cartesian)
           *res_map*
               Map resolution               
        """

        origin = densMap.origin
        apix = densMap.apix
        box_size = densMap.box_size()
        x_pos = int(round((atom.x-origin[0])/apix,0))
        y_pos = int(round((atom.y-origin[1])/apix,0))
        z_pos = int(round((atom.z-origin[2])/apix,0))
        #print atom.get_res_no(), x_pos,y_pos,z_pos
        if((densMap.x_size() > x_pos >= 0) and (densMap.y_size() > y_pos >= 0) and (densMap.z_size() > z_pos >= 0)):
            #search all points withing 1.5sigma
            list_points = gridtree.query_ball_point([atom.x,atom.y,atom.z], 1.5*max(sim_sigma_coeff*res_map,1.0))
            return list_points
        else:
            print("Warning, atom out of map box")
            return []
    def model_tree(self,list_coord1,distpot=6.0,list_coord2=None):
        """
        Returns 
        """
        try: 
            from scipy.spatial import cKDTree
            coordtree = cKDTree(list_coord1)
            if list_coord2 != None: 
                coordtree1 = cKDTree(list_coord2)
        except ImportError:
            from scipy.spatial import KDTree
            coordtree = KDTree(list_coord1)
            if list_coord2 != None: coordtree1 = KDTree(list_coord2)
        if list_coord2 != None: 
            neigh_points = coordtree.query_ball_tree(coordtree1,distpot)
            # use count_neighbors if the corresponding indices are not required
        else: 
            neigh_points = coordtree.query_ball_tree(coordtree,distpot)
        return neigh_points
        
    
    def get_coordinates(self,structure_instance):
        """
        
        Returns flat indices of the pixels occupied by each residue in a chain.
        
        Arguments:
        
           *structure_instance*
               Structure instance of the model.
        """
        dict_res_CA = {}
        dict_res_indices = {}
        dict_chain_indices = {}
        dict_chain_CA = {}    
        currentChain = structure_instance.atomList[0].chain
        for x in structure_instance.atomList:
            if not x.chain == currentChain:
                try: dict_chain_indices[x.model][currentChain] = dict_res_indices.copy()
                except KeyError: 
                    dict_chain_indices[x.model] = {}
                    dict_chain_indices[x.model][currentChain] = dict_res_indices.copy()
                try: dict_chain_CA[x.model][currentChain] = dict_res_CA.copy()
                except KeyError:
                    dict_chain_CA[x.model] = {}
                    dict_chain_CA[x.model][currentChain] = dict_res_CA.copy()
                currentChain = x.chain
                dict_res_indices = {}
                dict_res_CA = {}
            cur_chain = x.chain
            cur_res = x.get_res_no()
            #save residue coords
            if x.atom_name == 'CA': #x.fullid
                #CA coordinates
                dict_res_CA[cur_res] = [x.x, x.y, x.z]
            try: dict_res_indices[cur_res].append([x.x, x.y, x.z])
            except KeyError: dict_res_indices[cur_res] = [[x.x, x.y, x.z]]
        if not dict_chain_indices.has_key(currentChain):
            #
            try: dict_chain_CA[x.model][currentChain] = dict_res_CA.copy()
            except KeyError:
                dict_chain_CA[x.model] = {}
                dict_chain_CA[x.model][currentChain] = dict_res_CA.copy()
            try: dict_chain_indices[x.model][currentChain] = dict_res_indices.copy()
            except KeyError: 
                dict_chain_indices[x.model] = {}
                dict_chain_indices[x.model][currentChain] = dict_res_indices.copy()
            #dict_chain_indices[currentChain] = dict_res_indices.copy()
        return dict_chain_indices, dict_chain_CA   
    # get nearest grid indices for a residue atoms
    def get_indices(self,structure_instance,emmap,res_map,sim_sigma_coeff=0.187):
        """
        
        Returns flat indices of the pixels occupied by each residue in a chain.
        
        Arguments:
        
           *structure_instance*
               Structure instance of the model.
           *emmap*
               Map instance the model is to be placed on.
           *res_map*
               Resolution of the map
           *sim_sigma_coeff*
               Sigma factor used for blurring
        """

        dict_res_indices = {}
        dict_res_dist = {}
        dict_chain_res = {}
        dict_chain_indices = {}
        gridtree = self.maptree(emmap)
        points = []
        #get chain details
        currentChain = structure_instance.atomList[0].chain
        for x in structure_instance.atomList:
            if not x.chain == currentChain:
                #uniquify lists
                for el in dict_res_indices:
                    tmplist = dict_res_indices[el][:]
                    setlist = set(tmplist)
                    dict_res_indices[el] = list(setlist)
                dict_chain_indices[currentChain] = dict_res_indices.copy()
                currentChain = x.chain
                dict_res_indices = {}
            cur_chain = x.chain
            cur_res = x.get_res_no()
            #save residue numbers in order
            try: 
              if not cur_res in dict_chain_res[currentChain]: dict_chain_res[currentChain].append(cur_res)
            except KeyError: dict_chain_res[currentChain] = [cur_res]
            #return indices covered by gaussian blur
            if not res_map is None: points = self.mapGridPositions(emmap,x,gridtree,res_map,sim_sigma_coeff)
            #return indices covered by vdW radii
            else: points = self.mapGridPositions_vdw(emmap,x,gridtree)
            if len(points) == 0:
                dict_res_indices[cur_res] = []
                continue
            if x.atom_name == 'CA': #x.fullid
                #CA coordinates
                dict_res_dist[cur_res] = [x.x, x.y, x.z]
            # get points occupied by the residue
            if dict_res_indices.has_key(cur_res):
                dict_res_indices[cur_res].extend(points)
            else: dict_res_indices[cur_res] = points
        if not dict_chain_indices.has_key(currentChain):
            #uniquify lists
            for el in dict_res_indices:
                tmplist = dict_res_indices[el][:]
                setlist = set(tmplist)
                dict_res_indices[el] = list(setlist)
            dict_chain_indices[currentChain] = dict_res_indices.copy()
        return dict_chain_indices, dict_chain_res,dict_res_dist


    # get nearest grid indices for a residue atoms
    def _get_indices1(self,structure_instance,emmap,res_map,sim_sigma_coeff=0.187):
        """
        
        Returns flat indices of the pixels occupied by each residue in a chain.
        
        Arguments:
        
           *structure_instance*
               Structure instance of the model.
           *emmap*
               Map instance the model is to be placed on.
           *res_map*
               Resolution of the map
           *sim_sigma_coeff*
               Sigma factor used for blurring
        """

        dict_res_indices = {}
        dict_res_dist = {}
        dict_chain_indices = {}
        gridtree = self.maptree(emmap)
        points = []
        #get chain details
        currentChain = structure_instance.atomList[0].chain
        for x in structure_instance.atomList:
            if not x.chain == currentChain:
                dict_chain_indices[currentChain] = dict_res_indices.copy()
                currentChain = x.chain
                dict_res_indices = {}    
            cur_chain = x.chain
            cur_res = x.get_res_no()
            points = self.mapGridPositions(emmap,x,gridtree,res_map,sim_sigma_coeff)
            if len(points) == 0:
                dict_res_indices[int(cur_res)] = []
                continue
            if x.atom_name == 'CA': #x.fullid
                #CA coordinates
                dict_res_dist[cur_res] = [x.x, x.y, x.z]
            if dict_res_indices.has_key(int(cur_res)):
                dict_res_indices[int(cur_res)].extend(points)
            else: dict_res_indices[int(cur_res)] = points
        #uniquify lists
        for el in dict_res_indices:
            tmplist = dict_res_indices[el][:]
            setlist = set(tmplist)
            dict_res_indices[el] = list(setlist)
        return dict_res_indices, dict_res_dist

        
#this two can be merged and be a unique function that return either the density or 1
#added by PAP
    def make_atom_overlay_map(self, densMap, prot):
        
        """
        
        Returns a Map instance with atom masses superposed on it.
        
        Arguments:
        
           *densMap*
               an empty (all densities zero) Map instance to superpose the atoms onto.
           *prot*
               a Structure instance.
               
        """
        densMap = densMap.copy()
        for atom in prot.atomList:
            #print atom.atom_name
            #print vdw_radii
            pos = self.mapGridPosition(densMap, atom)
	    #print pos
            if pos:
                densMap.fullMap[pos[2]][pos[1]][pos[0]] += pos[3]
        return densMap

    def make_atom_overlay_mapB(self, densMap, prot):
        
        """
        
        Returns a Map instance with atom masses superposed on it.
        
        Arguments:
        
           *densMap*
               an empty (all densities zero) Map instance to superpose the atoms onto.
           *prot*
               a Structure instance.
               
        """
        densMap = densMap.copy()
        for atom in prot.atomList:
            pos = self.mapGridPosition(densMap, atom)
            print(pos)
            if pos:
                densMap.fullMap[pos[2]][pos[1]][pos[0]] += pos[3]
        return densMap
    #ADDED BY PAP
    def make_atom_overlay_map1(self, densMap, prot):
        
        """
        
        Returns a Map instance with atom locations recorded on the nearest voxel with a value of 1.
        
        Arguments:
           
           *densMap*
               an empty (all densities zero) Map instance to superpose the atoms onto.
           *prot*
               a Structure instance.
               
        """
        densMap = densMap.copy()
        densMap.fullMap = densMap.fullMap * 0
        for atom in prot.atomList:
            #print 'Atom name : ',atom.atom_name
            if (atom.atom_name == 'C' or atom.atom_name == 'N' or atom.atom_name == 'CA' or atom.atom_name == 'O' or atom.atom_name == 'CB'):
                pos = self.mapGridPosition(densMap, atom)
                #print 'overlay index', pos
                if pos:
                    densMap.fullMap[pos[2]][pos[1]][pos[0]] = 1
        return densMap
    
    def make_model_grid(self,prot,spacing,ca_only=False,densMap=False):
        if not densMap:
            densMap = self.protMap(prot, spacing)
            print("WARNING: Use StructureBlurrer.gaussian_blur_box() to blured a map with a user defined defined cubic box")
            #from here till newMap.fullMap*=0 are few line of code that create an empty map with the new A/px of 1
            #this replace the make_clash_map(apix) function. they do the job but they need to be replaced with something more rigorous
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        #get flat indices of grid
        nz,ny,nx = densMap.fullMap.shape
        zg,yg,xg = mgrid[0:nz,0:ny,0:nx]
        indi = zip(xg.ravel(), yg.ravel(), zg.ravel())
        gridtree = self.maptree(densMap)
        
        currentChain = prot.atomList[0].chain
        for x in prot.atomList:
            if not x.chain == currentChain:
                currentChain = x.chain
            cur_chain = x.chain
            cur_res = x.get_res_no()
            
            #CA only
            if ca_only and not x.atom_name == 'CA': #x.fullid
                #CA coordinates
                continue
            
            points,gridpoint = self.mapGridPositions_vdw(densMap,x,gridtree)
            if len(points) == 0:
                continue
            
            grid_indices = [] 
            for p in points:
                grid_indices.append(indi[p])
            x.grid_indices = grid_indices[:]
            
            array_indices = array(grid_indices)
            ind_arrxyz = transpose(array_indices)
            ind_arrzyx = (ind_arrxyz[2],ind_arrxyz[1],ind_arrxyz[0])
            
            densMap.fullMap[ind_arrzyx] = 1.0
        return densMap

        
        

    def gaussian_blur(self, prot, resolution, densMap=False, sigma_coeff=0.356, normalise=True,filename="None"):
        
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in reciprocal space.

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
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

            *filename*
                output name of the map file.
                
        """
        #densMap= your map if you want to compare prot blurred with an exisiting map.
        #Daven always use that so that it blurred based on the experiment box
        if not densMap:
            densMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print("WARNING: Use StructureBlurrer.gaussian_blur_box() to blured a map with a user defined defined cubic box")
            #from here till newMap.fullMap*=0 are few line of code that create an empty map with the new A/px of 1
            #this replace the make_clash_map(apix) function. they do the job but they need to be replaced with something more rigorous
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        ##newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        ##newMap.fullMap *= 0
        newMap = densMap.copy()
        newMap.fullMap = zeros((z_s, y_s, x_s))
        newMap.apix = (densMap.apix*densMap.x_size())/x_s
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        fou_map = fourier_gaussian(fftn(newMap.fullMap), sigma)
        newMap.fullMap = real(ifftn(fou_map))
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        newMap.filename=filename
        newMap.update_header
        return newMap

    #add IF
    def gaussian_blur_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, sigma_coeff=0.356, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in reciprocal space.
    
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
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

            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        ##newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        ##newMap.fullMap *= 0
        newMap = densMap.copy()
        newMap.fullMap = zeros((z_s, y_s, x_s))
        newMap.apix = (densMap.apix*densMap.x_size())/x_s
        sigma = sigma_coeff*resolution
        newMap = self.make_atom_overlay_map(newMap, prot)
        fou_map = fourier_gaussian(fftn(newMap.fullMap), sigma)
        newMap.fullMap = real(ifftn(fou_map))
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap
    
 
    def hard_sphere(self,prot,resolution, densMap=False, normalise=True,filename="None"):
             
        """
        
        Returns a Map instance based on a Hard Sphere model of a protein.
        Usefull for rigid fitting (Topf et al, 2008)

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
            *filename*
                output name of the map file.
                
        """
        gridpx=min(resolution/4., 3.5)
        if not densMap:
            densMap = self.protMap(prot, min(resolution/4., 3.5), resolution)
            print("WARNING: Use StructureBlurrer.hard_sphere() to create a map with a user defined defined cubic box")
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0

        newMap = self.make_atom_overlay_mapB(newMap, prot)
        #new_map.fullMap = laplace(self.fullMap)
        print(gridpx)
        #newMap.fullMap = uniform_filter(newMap.fullMap,size=gridpx,mode='constant',cval=0.0)
        newMap = newMap.resample_by_box_size(densMap.box_size())
        #newMap.fullMap=newMap
        #newMap = newMap.resample_by_box_size(gridpx)
        if normalise:
            newMap = newMap.normalise()
        return newMap
    #add IF
    def hard_sphere_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Hard Sphere model of a protein.
        Usefull for rigid fitting (Topf et al, 2008)
            
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        newMap.fullMap *= 0

        newMap = self.make_atom_overlay_map(newMap, prot)
        #newMap.fullMap=newMap
        #newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap

    #add IF
    def gaussian_blur_real_space(self, prot, resolution, densMap=False, sigma_coeff=0.356, normalise=True,filename="None"):
        
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in real space
        

        Arguments:

            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *densMap*
                False to build a Map with dimensions based on the protein, or a Map instance to be used as a template.
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

           *filename*
                output name of the map file.
                
        """
        if not densMap:
            #>
            densMap = self.protMap(prot, max(1.0,min(resolution/4., 3.5)), resolution)
            print("WARNING: Use StructureBlurrer.gaussian_blur_real_space_box() to blured a map with a user defined defined cubic box")
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        ##newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        ##newMap.fullMap *= 0
        newMap = densMap.copy()
        newMap.fullMap = zeros((z_s, y_s, x_s))
        newMap.apix = (densMap.apix*densMap.x_size())/x_s
        sigma = max(sigma_coeff*resolution,1.0)
        newMap = self.make_atom_overlay_map(newMap, prot)
        gauss_map = gaussian_filter(newMap.fullMap, sigma)
        newMap.fullMap = gauss_map
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap

    def gaussian_blur_real_space_box(self, prot, resolution, box_size_x, box_size_y, box_size_z, sigma_coeff=0.356, normalise=True,filename="None"):
        """
        
        Returns a Map instance based on a Gaussian blurring of a protein.
        The convolution of atomic structures is done in real space
           
        Arguments:
        
            *prot*
                the Structure instance to be blurred.
            *resolution*
                the resolution, in Angstroms, to blur the protein to.
            *box_size_x*
                 x dimension of map box in Angstroms.
            *box_size_y*
                y dimension of map box in Angstroms.
            *box_size_z* 
                z dimension of map box in Angstroms.
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

            *filename*
                output name of the map file.
                
        """
        densMap = self.protMapBox(prot, 1, resolution, box_size_x, box_size_y, box_size_z, filename)
        x_s = int(densMap.x_size()*densMap.apix)
        y_s = int(densMap.y_size()*densMap.apix)
        z_s = int(densMap.z_size()*densMap.apix)
        ##newMap = densMap.resample_by_box_size([z_s, y_s, x_s])
        ##newMap.fullMap *= 0
        newMap = densMap.copy()
        newMap.fullMap = zeros((z_s, y_s, x_s))
        newMap.apix = (densMap.apix*densMap.x_size())/x_s
        sigma = max(sigma_coeff*resolution,1.0)
        newMap = self.make_atom_overlay_map(newMap, prot)
        gauss_map = gaussian_filter(newMap.fullMap, sigma)
        newMap.fullMap = gauss_map
        newMap = newMap.resample_by_box_size(densMap.box_size())
        if normalise:
            newMap = newMap.normalise()
        return newMap


    
    #---BANDPASS FILTERING (NOT WORKING YET)--- add by DV# MAKE them PRIVITA _FUNCT
    #way of filtering the map using "Fourier-like" but it is too slow so abandon the idea. there are quiker and better way
    # Bsoft is a better way to go. http://lsbr.niams.nih.gov/bsoft/
    # not spend time on it.      
        

    def _bandpass_blur(self, atomList, densMap, lopass, lomin, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        pass
    

    def _bandpass_mask_gaussian(self, densMap, lopass, lopass_min, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        newMap = densMap.copy()#self.make_empty_map(densMap)
        centre = (array(newMap.box_size[:])-1)/2.0
        from time import time
        for z in range(newMap.box_size[2]):
            for y in range(newMap.box_size[1]):
                for x in range(newMap.box_size[0]):
                    t1 = time()
                    dist = sqrt((x-centre[0])**2 + (y-centre[1])**2 + (z-centre[2])**2)
                    t2 = time()
                    newMap[z][y][x] = self.bandpass_eq_gaussian(dist, lopass, lopass_min, lowid, hipass, hiwid)
                    t3 = time()
                    print(t2-t1, t3-t2)
        return newMap

    def _bandpass_eq_gaussian(self, dist, lopass, lopass_min, lowid, hipass, hiwid):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        lp_max = lopass+lowid
        hp_min = hipass-hiwid
        if dist <= lp_max:
            return lopass_min+(1-lopass_min)*exp(-0.5*((dist-lp_max)/lowid)**2)
        elif lp_max < dist <= hp_min:
            return 1.0
        else:
            return exp(-0.5*((dist-hp_min)/hiwid)**2)

    def _bandpass_test(self, lopass, lopass_min, lowid, hipass, hiwid, l_len):
        """
        
        WARNING: BANDPASS FILTERING (NOT WORKING YET)
        
        """
        from time import time
        start = time()
        a = zeros([l_len])
        for x in range(l_len):
            a[x] = self.bandpass_eq_gaussian(x, lopass, lopass_min, lowid, hipass, hiwid)
        end = time()
        print(end-start)
        return a
    
