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

from numpy import sqrt,matrix,random
from math import cos,sin,pi,acos,asin, atan2, radians, sqrt as msqrt, pow as mpow

class Vector:
    """A class representing Cartesian 3-dimensonal vectors."""
    
    def __init__(self, x,y,z):
        """x, y, z = Cartesian co-ordinates of vector."""
        self.x = x
        self.y = y
        self.z = z
        
    def __repr__(self):
        return "(%.3f,%.3f,%.3f)" %(self.x, self.y, self.z)

    def __getitem__(self, index):
        l = [self.x, self.y, self.z]
        return l[index]

    def __iter__(self):
        l = [self.x, self.y, self.z]
        return l.__iter__()

    def copy(self):
        """
        Return:
            A copy of Vector instance
        """
        return Vector(self.x, self.y, self.z)
    
    def mod(self):
        """
        Return:
            The modulus (length) of the vector.
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def reverse(self):
        """
        Flip the direction of a Vector instance.
        
        Return:
          A Vector instance  
        """
        return Vector(-self.x,-self.y,-self.z)
    
    def arg(self, vector):
        """Return the argument (angle) between this and another vector.RAD"""
        top  = self.dot(vector)
        bottom = self.mod()*vector.mod()
        #print top, bottom
        #print 'top/bottom, top, bottom ', top/bottom, top, bottom
        if abs(top-bottom) < 0.00001:
            return 0.0
        else:
            #print 'top/bottom, top, bottom ', top/bottom, top, bottom
            #return acos(top/bottom)
            return acos(min(max(top/bottom,-1.0),1.0))
    
    def times(self, factor):
        """
        Multiplies a Vector instance by a scalar factor.
        
        Return:
          A Vector instance
        """
        return Vector(factor*self.x, factor*self.y, factor*self.z)
    
    def dot(self, vector):
        """
        Return:
            The dot product of this and another vector specified as input parameter. 
        """
        return vector.x * self.x + vector.y * self.y + vector.z * self.z
    
    def cross(self, vector):
        """
        Return:
            A Vector instance of the cross product of this and another vector specified as input parameter
        """
        newX = self.y*vector.z - self.z*vector.y
        newY = self.z*vector.x - self.x*vector.z
        newZ = self.x*vector.y - self.y*vector.x
        return Vector(newX, newY, newZ)
    
    def __sub__(self, vector):
        """Return a Vector instance of the subtraction of a vector from this one."""
        newX = self.x - vector.x
        newY = self.y - vector.y
        newZ = self.z - vector.z
        return Vector(newX, newY, newZ)

    def dist(self, vector):
        """
        Return:
            The distance between this and another vector specified as input parameter. 
        """
        return (self-vector).mod()
    
    def __add__(self, vector):
        """Return a Vector instance of the addition of a vector from this one."""
        newX = self.x + vector.x
        newY = self.y + vector.y
        newZ = self.z + vector.z
        return Vector(newX, newY, newZ)

    def __mul__(self, prod):
        newX = self.x*prod
        newY = self.y*prod
        newZ = self.z*prod
        return Vector(newX, newY, newZ)

    def __div__(self, divisor):
        newX = self.x/float(divisor)
        newY = self.y/float(divisor)
        newZ = self.z/float(divisor)
        return Vector(newX, newY, newZ)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def translate(self, x, y, z):
        """
        Translate a Vector instance.
        
        Arguments:
            *x, y, z*
                distance in Angstroms in respective Cartesian directions to translate vector.
                
        Return:
            Vector instance.    
        """
        newX = self.x + x
        newY = self.y + y
        newZ = self.z + z
        return Vector(newX, newY, newZ)

    def matrix_transform(self, rot_mat):
        """
        Transform the vector using a transformation matrix.
        
        Arguments:
            *rot_mat*
                a 3x3 Python matrix instance.
        Return:
            A vector instance
        """
        vec_mat = matrix([[self.x],[self.y],[self.z]])
        new_pos = rot_mat*vec_mat
        x = float(new_pos[0])
        y = float(new_pos[1])
        z = float(new_pos[2])
        return Vector(x,y,z)

    def to_atom(self):
        """
        Create an Atom instance based on Vector instance.
        
        Return:
            Atom instance
        """
        from ProtRep_Biopy import BioPyAtom
        template = 'ATOM      1  C   NOR A   1      23.161  39.732 -25.038  1.00 10.00             C'
        a = BioPyAtom(template)
        a.x = self.x
        a.y = self.y
        a.z = self.z
        return a
    
    def unit(self):
        """
        Return:
            Vector instance of a unit vector.
        """
        mod = self.mod()
        if mod==0:
            return Vector(0,0,0)
        return Vector(self.x/mod, self.y/mod, self.z/mod)

###########################################################################
###########################################################################
###########################################################################
#### def out of the class . 
#### better have them separate as these definition are an adaptation of 
#### Transformations Python Module from Christoph Gohlke
#### http://www.lfd.uci.edu/~gohlke/
#### http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
###########################################################################
###########################################################################
###########################################################################
   
def random_vector(min_v, max_v):
    """
    Generate a random vector.
    The values for the vector component x, y, and z are randomly sampled between minimum and maximum values specified.
    
    Argument:
        *min_v, max_v*
            minimum and maximum value
    Return:
        A Vector instance.
    """
    
    x = random.uniform(min_v, max_v)
    y = random.uniform(min_v, max_v)
    z = random.uniform(min_v, max_v)
    return Vector(x,y,z)

def axis_angle_to_matrix(x, y, z, turn, rad=False):
    """
    Converts the axis angle rotation to a matrix form.
    
    Arguments:
       *x, y, z*
           axis of rotation (does not need to be normalised).
       *turn*
           angle of rotation, in radians if rad=True, else in degrees.
    Return:
        A 3X3 transformation matrix.   
    """
    if not rad:
        turn = turn*pi/180
    c_a = cos(turn)
    s_a = sin(turn)
    v = Vector(x,y,z).unit()
    x = v.x
    y = v.y
    z = v.z

    rot_mat = matrix([[x**2+(1-x**2)*c_a, x*y*(1-c_a)-z*s_a, x*z*(1-c_a)+y*s_a],
                          [ x*y*(1-c_a)+z*s_a, y**2+(1-y**2)*c_a, y*z*(1-c_a)-x*s_a],
                          [x*z*(1-c_a)-y*s_a, y*z*(1-c_a)+x*s_a, z**2+(1-z**2)*c_a]])
    return rot_mat

def euler_to_matrix(x_turn, y_turn, z_turn, rad=False):
    """
    Converts an euler rotation to a matrix form.
    
    Arguments:
       *x_turn, y_turn, z_turn*
           rotation angles around respective axis, in radians if rad=True, else in degrees.
    Return:
        A 3X3 transformation matrix.   
    """
    if not rad:
        x_turn = x_turn*pi/180
        y_turn = y_turn*pi/180
        z_turn = z_turn*pi/180
    x_mat = axis_angle_to_matrix(0,0,1,x_turn, rad=True)
    y_mat = axis_angle_to_matrix(0,1,0,y_turn, rad=True)
    z_mat = axis_angle_to_matrix(1,0,0,z_turn, rad=True)
    return x_mat*y_mat*z_mat

def axis_angle_to_euler(x,y,z, turn, rad=False):
    """
    Converts the axis angle rotation to an Euler form.
    
    Arguments:
       *x, y, z*
           axis of rotation (does not need to be normalised).
       *turn*
           angle of rotation, in radians if rad=True, else in degrees.
    Returns:
        A 3-tuple (x,y,z) containing the Euler angles. .
    """
    if not rad:
        turn = turn*pi/180
    z_rot = atan2(y*sin(turn)-x*z*(1-cos(turn)), 1-(y**2+z**2)*(1-cos(turn)))
    x_rot = asin(x*y*(1-cos(turn))+z*sin(turn))
    y_rot = atan2(x*sin(turn)-y*z*(1-cos(turn)), 1-(x**2+z**2)*(1-cos(turn)))
    return (x_rot, y_rot, z_rot)

# -- Vector methods for torsion angle geometry -- #



def torsion(a, b, c):
    """
    Find the torsion angle between planes ab and bc.
    
    Arguments:
        
        *a,b,c*
            Vector instances.
    
    Returns:
        The torsion angle in radians
    """
    n1 = a.cross(b)
    n2 = b.cross(c)
    return n1.arg(n2)

def calcMtrx(arr):
    """
    Calculate 3 x 4 transformation matrix from Euler angles and offset.    
    Arguments:        
        *arr*
            [psi,theta,phi,offsetx,offsety,offsetz].    
    Returns:
        3 x 4 transformation matrix
    """
    # Taken from http://mathworld.wolfram.com/EulerAngles.html 
    # Order psi/theta/phi often reversed
    Cpsi = cos(radians(float(arr[0])))
    Spsi = sin(radians(float(arr[0])))
    Ctheta = cos(radians(float(arr[1])))
    Stheta = sin(radians(float(arr[1])))
    Cphi = cos(radians(float(arr[2])))
    Sphi = sin(radians(float(arr[2])))

    res = [[0 for row in range(4)] for col in range(4)]

    res[0][0] = Cpsi * Cphi - Ctheta * Sphi * Spsi
    res[0][1] = Cpsi * Sphi + Ctheta * Cphi * Spsi
    res[0][2] = Spsi * Stheta
    res[1][0] = -Spsi * Cphi - Ctheta * Sphi * Cpsi
    res[1][1] = -Spsi * Sphi + Ctheta * Cphi * Cpsi
    res[1][2] = Cpsi * Stheta
    res[2][0] = Stheta * Sphi
    res[2][1] = -Stheta * Cphi
    res[2][2] = Ctheta
    '''
    # y convention
    res[0][0] = -Spsi * Sphi + Ctheta * Cphi * Cpsi
    res[0][1] = Spsi * Cphi + Ctheta * Sphi * Cpsi
    res[0][2] = -Cpsi * Stheta
    res[1][0] = -Cpsi * Sphi - Ctheta * Cphi * Spsi
    res[1][1] = Cpsi * Cphi - Ctheta * Sphi * Spsi
    res[1][2] = Spsi * Stheta
    res[2][0] = Stheta * Cphi
    res[2][1] = Stheta * Sphi
    res[2][2] = Ctheta
    '''

    res[0][3] = float(arr[3])
    res[1][3] = float(arr[4])
    res[2][3] = float(arr[5])
    return res                         

def _rotmat_to_axisangle(mat):
    """
    Convert rotation matrix to axisangle.
    Arguments:
        *mat*
            Rotation matrix.
    Returns:
        axis-angle
    """
    trace_mat=float(mat[0,0])+float(mat[1,1])+float(mat[2,2])

    if trace_mat > 3.0:
        print '_rotmat_to_axisangle trace_mat to large:', trace_mat
        trace_mat =  3.0
    if trace_mat < -1.0:
        print '_rotmat_to_axisangle trace_mat to small:', trace_mat
        trace_mat = -1.0

    theta1=acos((trace_mat-1)/2.0)
    if theta1 == 0.0: return (0.0, 1.0,1.0,1.0)
    k=1/(2.0*sin(theta1))
    # vector axis:
    n1=k*(mat[2,1]-mat[1,2])
    n2=k*(mat[0,2]-mat[2,0])
    n3=k*(mat[1,0]-mat[0,1])
    return (theta1,n1,n2,n3)

def cps(mat_1,mat_2):
    """
    Find rotation and translation difference between two transformations.
    Arguments:
        *mat_1,mat_2*
            Transformation matrices.
    Returns:
        The translation and rotation differences
    """
    mat1 = matrix(mat_1)
    mat2 = matrix(mat_2)

    # mat to euler (angular distance not accurate!)
    #t1,p1,s1 = _rotmat_to_euler(mat1)
    #t2,p2,s2 = _rotmat_to_euler(mat2)
    #ang_magnitude = msqrt(mpow(t2-t1,2)+mpow(p2-p1,2)+mpow(s2-s1,2))
    matR1 = matrix([mat_1[0][:-1],mat_1[1][:-1],mat_1[2][:-1]])
    matR2 = matrix([mat_2[0][:-1],mat_2[1][:-1],mat_2[2][:-1]])
    matR = matR2*matR1.transpose()
    ang_magnitude,xV,yV,zV = _rotmat_to_axisangle(matR)
    ang_magnitude = ang_magnitude*(180.0/pi)
    shift = msqrt(mpow(mat2[0,3]-mat1[0,3],2)+mpow(mat2[1,3]-mat1[1,3],2)+mpow(mat2[2,3]-mat1[2,3],2))
    #print (mpow(mat2[0,3]-mat1[0,3],2)+mpow(mat2[1,3]-mat1[1,3],2)+mpow(mat2[2,3]-mat1[2,3],2)), ang_magnitude
    #acps_score = (pi/360)*(mpow(mat2[0,3]-mat1[0,3],2)+mpow(mat2[1,3]-mat1[1,3],2)+mpow(mat2[2,3]-mat1[2,3],2))*abs(ang_magnitude)
    return shift, ang_magnitude


def altTorsion(a,b,c):
    """
    An alternate and better way to find the torsion angle between planes ab and bc.
    
    Arguments:
        *a,b,c*
            Vector instances.
    Return:
        The torsion angle (radians)
    """
    A = a.dot(b.cross(c))*b.mod()
    B = (a.cross(b)).dot(b.cross(c))
    return atan2(A,B)
#PAP addition
def random_vector2(ul_list):
    x = random.uniform(ul_list[0],ul_list[1])
    y = random.uniform(ul_list[0],ul_list[1])
    z = random.uniform(ul_list[0],ul_list[1])
    return Vector(x,y,z)

def align_2seqs(seq1,seq2):
  try:
    #NW from biopython
    import Bio.pairwise2
    aln1 = Bio.pairwise2.align.globalms(''.join(seq1),''.join(seq2),1,-1,-0.5,-0.1)[0]
    align_seq1,align_seq2 = list(aln1[0]),list(aln1[1])
  except:
    print 'ImportError: Bio.pairwise2'
    return
  return align_seq1,align_seq2