##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright 2010-2014 TEMPy Inventors and Birkbeck College University of London.
#                          The TEMPy Inventors are: Maya Topf, Daven Vasishtan, 
#                           Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph,
#                          Harpal Sahota
# 
# 
#     TEMPy is available under Public Licence.
#     
#     Please cite your use of TEMPy in published work:
#     
#     Vasishtan D, Topf M. (2011) J Struct Biol 174:333-343. Scoring functions for cryoEM density fitting.
#     Pandurangan AP, Vasishtan D, Topf M. (2015) Structure 23:2365-2376. GAMMA-TEMPy: Simultaneous fitting of components in 3D-EM Maps of their assembly using genetic algorithm.
#===============================================================================

from TEMPy.MapParser import *
from TEMPy.ProtRep_Biopy import BioPy_Structure
from TEMPy.Vector import *
from numpy import random, nonzero, arange, ones, zeros, repeat, sort, power
from time import time
import numpy as np
import math
def VQ(D, n, epochs, alpha0=0.5, lam0=False):
    """ 
    
    Function to clusters a set of vectors (D) into a number (n) of codebook vectors

        Arguments:
            *n*
                number of VQ points to generate
            *epochs*
                number of iterations to run the algorithm
            *alpha0* 
                NOTE: Ask Daven regarding this parameter.
            *lam0* 
                NOTE: Ask Daven regarding this parameter.
    
    """

    if not lam0:
        lam0 = n/2.

    dlen, dim = D.shape

    neurons = (np.random.rand(n, dim-1)-0.5)*0.2
    train_len = epochs*dlen
    #random.seed(int(100*time()))

    den_sum = sum(D[:,3])

    den_culm = []
    den_culm.append(D[0,3])
    for i in range(1,dlen):
        den_culm.append(den_culm[i-1]+D[i,3])
    den_culm = np.matrix(den_culm)
    
    den_rand = den_sum*np.random.rand(train_len,1)
    sample_inds = []
    for j in range(train_len):
        #if j%100000 == 0:
        #    print str(j)+' j'
        min_index = den_rand[j]<den_culm
        sample_inds.append(np.nonzero(min_index)[1][0,0])
    lam = lam0*(0.01/lam0)**(np.arange(train_len)/float(train_len))
    alpha = alpha0*(0.005/alpha0)**(np.arange(train_len)/float(train_len))

    for i in range(train_len):
        #if i%5000 == 0:
        #    print i
        x = D[sample_inds[i],:3]
        x = np.array([int(y) for y in x])
        known = np.array([not math.isnan(y) for y in x])
        X = np.repeat([x], [n], axis=0)#x[ones((n,1), dtype='i'), known]
        X = np.array(X[:,known], dtype='i')

        Dx = np.matrix(neurons[:, known] - X)
        inds = sorted(enumerate(np.power(Dx,2)*np.matrix(known).T), key=lambda t: t[1])
        inds = np.array([a[0] for a in inds])
        ranking = np.zeros(len(inds))
        ranking[inds] = range(n)

        h = map(math.exp, -ranking/lam[i])
        H = np.repeat([h], [len(known)], axis=0).T

        neurons = neurons + (alpha[i]*H)*(X-neurons[:, known])
    return [Vector(x[0], x[1], x[2]) for x in neurons]

def map_points(emmap, threshold):
    """ 
    
    Function to return all map grid coordinated whose density value is greater than a given cutoff value.

        Arguments:
            *emmap*
                Instance of a map object.
            *threshold*
                Density value used as a cutoff.
    
    """

    emvec = np.nonzero(emmap.fullMap>threshold)
    apix = emmap.apix
    x_o = emmap.x_origin()
    y_o = emmap.y_origin()
    z_o = emmap.z_origin() 
    emv = []
    for v in range(len(emvec[0])):
        x = emvec[2][v]*apix+x_o
        y = emvec[1][v]*apix+y_o
        z = emvec[0][v]*apix+z_o
        dens = emmap[emvec[0][v], emvec[1][v], emvec[2][v]]
        emv.append([x,y,z,dens])
    return np.array(emv)

def write_to_pdb(vq, output_file=None):
    """ 
    
    Function to write the VQ point in PDB ATOM record format.

        Arguments:
            *vq*
                A list of vector object containing the VQ points.
            *output_file*
                Name of the output file.
    
    """

    f = open(output_file+'_vq.pdb', 'w')
    header='''VQ POINTS GENERATED WITH TEMPY'''
    record_name = 'ATOM'
    serial = 1
    atom_name = 'CA'
    alt_loc = ' '
    chain = 'A'
    res_no = 1
    res = 'GLY'
    icode = ''
    occ = 1.0
    temp_fac = 1.0
    elem = 'C'
    charge = ''
    for v in vq:
	line = ''
        line += record_name.ljust(6) 
        line += str(serial).rjust(5)+' '
        line += atom_name.center(4)
        line += alt_loc.ljust(1)
        line += res.ljust(3)+' '
        line += chain.ljust(1)
        line += str(res_no).rjust(4)
        line += str(icode).ljust(1)+'   '
        x = '%.3f' % v.x
        y = '%.3f' % v.y
        z = '%.3f' % v.z
        line += x.rjust(8)
        line += y.rjust(8)
        line += z.rjust(8)
        occ = '%.2f'% float(occ)
        temp_fac = '%.2f'% float(temp_fac)
        line += occ.rjust(6)
        line += temp_fac.rjust(6)+'          '
        line += elem.strip().rjust(2)
        line += charge.strip().ljust(2)
        line = line + '\n'
        f.write(line)
	serial += 1
	res_no += 1
    f.close()	
    #atomList = [x.to_atom() for x in vq]
    ##s = Structure(atomList)
    #s = BioPy_Structure(atomList)
    #if output_file:
    #    s.write_to_PDB(output_file)
    #return s

def get_VQ_points(emmap, threshold, noOfPoints, epochs, output_file=None, lap_fil=True):
    """
       Function to generate VQ point given a density map

       Arguments:

           *emmap* 
               Map (to be clustered) instance.
           *threshold* 
               voxels with density above this value are used in the VQ run.
           *noOfPoints* 
               num of VQ points to output.
           *epochs*
               num of iterations to run the algorithm
           *output_file* 
               file to output to. In PDB format
           *lap_fil* 
               True if you want to Laplacian filter the map first, False otherwise. Note that filtering the map 
               will change the density values of the map, which is relevant for the threshold parameter.

    """
    #ASK SHIHUA WHAT VALUE TO USE HERE for epochs
    if lap_fil:
        emmap = emmap.laplace_filtered()
    	emmap = emmap.normalise()
    D = map_points(emmap, threshold)
    vq = VQ(D, noOfPoints, epochs)
    if output_file:
        write_to_pdb(vq, output_file)
    #else:
    #    return vq
    return vq
