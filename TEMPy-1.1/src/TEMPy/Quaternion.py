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

from math import sin, tan, cos, acos, atan2, asin, sqrt, pi
from numpy import matrix
from random import randrange, uniform

class Quaternion:
    """ 
    
    A class representing quaternions.
    
    """

    
    def __init__(self, q_list):
        """
        
        A constructor to initialise the quaternion.
        
        
        Arguments:
            *q_list*
                A list containing data to load a quaternion. The format of the list is [w,x,y,z]
               
        """

	self.param = q_list
        
    def __repr__(self):
        """
        
        Return a string containing the details of the components of the quaternion.
        
        """

	return str(self.param)

    def copy(self):
        """
        
        Return an instance of the Quaternion object.
        
        """

        return Quaternion(self.param)
    
    def normalise(self, tolerance=0.00001):
        """
        
        Return a normalised quaternion vector.
        
        """

	mag1 = sum(n*n for n in self.param)
	if abs(mag1 - 1.0) > tolerance:
		mag2 = sqrt(mag1)
		q = list(n/mag2 for n in self.param)
	return Quaternion(q)

    def unit_quat(self):
        """
        
        Return a unit quaternion.
        
        """
	s = sum(n*n for n in self.param)
	mag = sqrt(s)
	return Quaternion(list(n/mag for n in self.param))

    def __mul__(self, other_q):
        """
        
        Return a quaternion object which the product of two quaternion.
 
        Arguments:
            *other_q*
                A list of type [w,x,y,z] to represent a quaternion vector.
        
        """

	q0, q1, q2, q3 = self.param
	r0, r1, r2, r3 = other_q.param
	w = r0*q0 - r1*q1 - r2*q2 - r3*q3
	x = r0*q1 + r1*q0 - r2*q3 + r3*q2
	y = r0*q2 + r1*q3 + r2*q0 - r3*q1
	z = r0*q3 - r1*q2 + r2*q1 + r3*q0
        return Quaternion([w,x,y,z])

    def multiply_3(self, obj1, obj2, obj3):
        """
        
        Return a quaternion object which the product of three quaternion.
 
        Arguments:
            *obj1*
                A list of type [w,x,y,z] to represent a quaternion vector.
            *obj2*
                A list of type [w,x,y,z] to represent a quaternion vector.
            *obj3*
                A list of type [w,x,y,z] to represent a quaternion vector.
        
        """

	q0, q1, q2, q3 = obj1.param
	r0, r1, r2, r3 = obj2.param
	w = r0*q0 - r1*q1 - r2*q2 - r3*q3
	x = r0*q1 + r1*q0 - r2*q3 + r3*q2
	y = r0*q2 + r1*q3 + r2*q0 - r3*q1
	z = r0*q3 - r1*q2 + r2*q1 + r3*q0

	q0, q1, q2, q3 = obj3.param
	r0, r1, r2, r3 = w, x, y, z
	w = r0*q0 - r1*q1 - r2*q2 - r3*q3
	x = r0*q1 + r1*q0 - r2*q3 + r3*q2
	y = r0*q2 + r1*q3 + r2*q0 - r3*q1
	z = r0*q3 - r1*q2 + r2*q1 + r3*q0

        return Quaternion([w,x,y,z])

    def conjuate(self, q):
        """
        
        Return a conjugate of a quaternion.
 
        Arguments:
            *q*
                A list of type [w,x,y,z] to represent a quaternion vector. NOTE: The argument seem to be not used. Will have to be removed and tested.
        
        """

	q0, q1, q2, q3 = self.param
	return Quaternion([q0, -q1, -q2, -q3])

    def mag(self):
        """
        
        Return the magnitude of the quaternion.
 
        
        """

	q0, q1, q2, q3 = self.param
	return sqrt(q0*q0+q1*q1+q2*q2+q3*q3)

    def to_rotation_matrix(self):
        """
        
        Convert the quaternion vector to a rotation matrix and returns the rotation matrix.
 
        
        """

	a, b, c, d = self.param
	rot_mat = matrix([[a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*a*c+2*b*d],
                          [2*b*c+2*a*d, a*a-b*b+c*c-d*d, 2*c*d-2*a*b],
                          [2*b*d-2*a*c, 2*a*b+2*c*d, a*a-b*b-c*c+d*d]])
	return rot_mat
