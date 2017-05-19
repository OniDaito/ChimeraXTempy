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

import numpy
from numpy import array, fromfile, flipud,isnan, transpose
import struct as binary
import string
from TEMPy.EMMap import Map

class MapParser:
    """
    A class to read various EM map file types into a Map object instance.
    """

    def __init__(self):
        ## mapping of numpy type to MRC mode
        self.numpy2mrc = {
            ## convert these to int8
            numpy.uint8: 0,
            numpy.bool: 0,
            numpy.bool_: 0,

            ## convert these to int16
            numpy.int16: 1,
            numpy.int8: 1,

            ## convert these to float32
            numpy.float32: 2,
            numpy.float64: 2,
            numpy.int32: 2,
            numpy.int: 2,

            ## convert these to complex64
            numpy.complex: 4,
            numpy.complex64: 4,
            numpy.complex128: 4,

            ## convert these to uint16
            numpy.uint16: 6,
        }

    @staticmethod
    def readMRCHeader(filename, endian = '<'):
        """
        Gets the header information from the MRC map file. 
        
        Argument
           *filename*
               input MRC map file name
           *endian*
               Endianness: Little or big
                
        Return:
           A string containing the MRC header information.

        """
        f = open(filename,'rb')
        fm_string = endian+(10*'l')+(6*'f')+(3*'l')+(3*'f')+(27*'l')+(3*'f')+(4*'c')+'lfl'
        header = list(binary.unpack(fm_string, f.read(224)))
        notes = f.read(800)
        notes = string.replace(notes, '\x00', '')
        header.append(notes)
        header = tuple(header)
        f.close()
        return header

    @staticmethod
    def get_endian(filename):
        """
        Read an MRC map file
           
        Arguments:
            *filename* 
                input MRC map file name.
        
        Return:
           Endianness: Little or big
        """
        h = MapParser.readMRCHeader(filename)
        if 0 <= h[3] <= 6:
            endian = '<'
        else:
            endian = '>'
        return endian

    @staticmethod
    def readMRC(filename):
        """
        Read an MRC map file
           
        Arguments:
            *filename* 
                input MRC map file name.
        
        Return:
           A Map instance containing the data read from MRC map file.
        """
        
        mrc2numpy = {
            0: numpy.uint8,
            1: numpy.int16,
            2: numpy.float32,
            #    3:  complex made of two int16.  No such thing in numpy
            #   however, we could manually build a complex array by reading two
            #   int16 arrays somehow.
            4: numpy.complex64,
            6: numpy.uint16,    # according to UCSF
        }

        endian = MapParser.get_endian(filename)
        header = MapParser.readMRCHeader(filename, endian)

        box_size = tuple(flipud(header[0:3]))
        origin = header[49:52] #ctrl UCSF

        # READ ORIGIN BASED ON MRC2000/CCP4 format
        nstart_index = header[4:7]
        apix = header[10]/header[0]
        nstart = (header[4]*float(apix),header[5]*float(apix),header[6]*float(apix))
        crs_index = header[16:19]
        if not (1 in (crs_index[0], crs_index[1], crs_index[2]) and 2 in (crs_index[0], crs_index[1], crs_index[2]) and 3 in (crs_index[0], crs_index[1], crs_index[2])):
                crs_index = (1,2,3)
        #print 'Axis order: ', crs_index
        #print 'Nstart', nstart_index[crs_index[0]-1],nstart_index[crs_index[1]-1],nstart_index[crs_index[2]-1]

	flag_orig = 0
        list_orig = [0.0, 0.0, 0.0]

        try:
                if header[52:56] == ('M','A','P',' '):
                        #print 'MAP flag found (MRC2000)'
                        origin = header[49:52]
                        #print 'Origin record: ', origin
                        if (isnan(origin[0]) or isnan(origin[1]) or isnan(origin[2])) or (origin[0] == 0.0 and origin[1] == 0.0 and origin[2] == 0.0):
                                origin = (0.0, 0.0, 0.0)
                                #print 'ORIGIN record empty, Checking NSTART records'
                                flag_orig = 1
                else:
                        flag_orig = 1
        except IndexError:
                origin = (0.0, 0.0, 0.0)
                pass

        if flag_orig == 1:
                if (nstart[0] == 0 and nstart[1] == 0 and nstart[2] == 0) or (isnan(nstart[0]) or isnan(nstart[1]) or isnan(nstart[2])):
                        #print 'NSTART records empty'
                        origin = (0.0, 0.0, 0.0)
                else:
                        list_orig[crs_index[0]-1] = nstart[0]
                        list_orig[crs_index[1]-1] = nstart[1]
                        list_orig[crs_index[2]-1] = nstart[2]
                        origin = (list_orig[0],list_orig[1],list_orig[2])

	'''
        if (nstart[0] == 0 and nstart[1] == 0 and nstart[2] == 0) or (isnan(nstart[0]) or isnan(nstart[1]) or isnan(nstart[2])):
                origin = (0.0, 0.0, 0.0)
                #print 'NSTART records empty, Checking ORIGIN records...'
                try:
                        if header[52:56] == ('M','A','P',' '):
                                #print 'MAP flag found (MRC2000)'
                                origin = header[49:52]
                                #print 'Origin record: ', origin
                                if (isnan(origin[0]) or isnan(origin[1]) or isnan(origin[2])) or (origin[0] == 0.0 and origin[1] == 0.0 and origin[2] == 0.0):
                                        origin = (0.0, 0.0, 0.0)
                except IndexError: pass
        else:
                list_orig[crs_index[0]-1] = nstart[0]
                list_orig[crs_index[1]-1] = nstart[1]
                list_orig[crs_index[2]-1] = nstart[2]
                origin = (list_orig[0],list_orig[1],list_orig[2])

        #print 'Map Origin is: ', origin
	'''

        map_size = header[0]*header[1]*header[2]
        f = open(filename,'rb')
        f.seek(1024)
        map_data = fromfile(f, dtype=mrc2numpy[header[3]], count=map_size)
        ### Swap bytes for endian
        if endian == '>':
            #print 'Byte order swapped!'
            map_data.byteswap(True)
	map_data=map_data.reshape(box_size)
        map_data=array(map_data, dtype='float64')

	### Check crs to xyz match
        if crs_index[0] != 1 or crs_index[1] != 2 or crs_index[2] != 3:
                #print 'Map axis permuted!!'
                #crs to xyz
                list_ind = [crs_index[0]-1,crs_index[1]-1,crs_index[2]-1]
                #xyz to crs
                index_new = (list_ind.index(0),list_ind.index(1),list_ind.index(2))
                #rearrange
                index_new1 = [2-index_new[2-a] for a in (0,1,2)]
                map_data=transpose(map_data,index_new1)

        f.close()
        return Map(map_data, origin, apix, filename, header=header)

    #BROKEN
    @staticmethod
    def _readXPLOR(filename, user_origin=None, user_box_size=None):
        """
        Read density map file in XPLOR format
        NOTE: broken.
        
        Argument:
           *filename*
               input XPLOR map file name.
        
        Return:
           A Map instance containing the data read from XPLOR map file.

        """
        f = open(filename, 'r')
        while(True):
            l = f.readline().split()
            #print l
            if(len(l) ==1 and l[0] == '0'):
                break
        new_map = []
        line = 1
        while(True):
            line = f.readline().split()
            for dens in line:
                new_map.append(float(dens))
            if len(new_map) >= new_map.box_size[0]*new_map.box_size[1]*new_map.box_size[2]:
                break
        new_map = array(new_map)
        new_map = new_map.reshape(new_map.box_size[2], new_map.box_size[1], new_map.box_size[0])
        f.close()
        return Map(new_map, new_map.origin, new_map.box_size, new_map.apix)

    @staticmethod
    def _readSitus(self,filename):
        """
        Read density map file in Situs format
        
        Arguments:
           *filename*
               input Situs map file name.
        
        Return:
            A Map instance containing the data read from Situs map file.

        """
        f = open(self.filename, 'r')
        first_line = f.readline().split()
        apix = float(first_line[0])
        origin = map(float, first_line[1:4])
        box_size = map(int, first_line[4:7])
        new_map = []
        line = 1
        while(True):
            line = f.readline().split()
            for dens in line:
                new_map.append(float(dens))
            if len(new_map) >= box_size[0]*box_size[1]*box_size[2]:
                break
        new_map = array(new_map)
        new_map = new_map.reshape(box_size[2], box_size[1], box_size[0])
        f.close()
        return Map(new_map, origin, box_size, apix)
