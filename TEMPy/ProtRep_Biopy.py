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

from numpy import ndarray,array,append,matrix
from math import pi
from random import randrange
import TEMPy.Vector as Vector
import sys
from TEMPy.Quaternion import *


# Useful global constants
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ12345676890abcdefghijklmnopqrstuvwxyz'
atomicMasses = {'H':1, 'C': 12, 'N':14, 'O':16, 'S':32}
vdw_radii = {'H':1.09, 'C': 1.7, 'N':1.55, 'O':1.52, 'S':1.8, 'P':1.8, 'Cl':1.75, 'CL':1.75, 'Cu':1.4, 'CU':1.4} # Taken from http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4
sequenceConsts = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I', 'MET':'M', 'PHE':'F', 'TRP':'W',\
                  'PRO':'P', 'SER':'S', 'THR':'T', 'CYS':'C', 'TYR':'Y', 'ASN':'N', 'GLN':'Q', 'ASP':'D',\
                  'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H'}
root2Pi = (2*pi)**0.5
aa_mass = {'A':71.0788,'R':156.1875,'N':114.1038,'D':115.0886,'C':103.1388,'E':129.1155,'Q':128.1307,'G':57.0519,'H':137.1411,'I':113.1594,'L':113.1594,'K':128.1741,'M':131.1926,'F':147.1766,'P':97.1167,'S':87.0782,'T':101.1051,'W':186.2132,'Y':163.1760,'V':99.1326}


class BioPyAtom:
    """
    
    A class representing an atom, as read from a PDB file using Biopython.
    
    """

    def __init__(self, atom):
        """Atom from BioPython"""
        if atom == []:
            return

        #http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
        #print "bioatom",atom#'bioatom <Atom O>'
        if atom.get_parent().get_id()[0][0] == "W" or atom.get_parent().id[0][0]=="H":
            self.record_name = "HETATM"
        else:
            self.record_name = "ATOM" # was pdbString[:6].strip() as "ATOM"
#             res.id[0] == "W" or res.id[0][0]=="H": #skip water and hetero residues
        self.serial = atom.get_serial_number()
        self.atom_name = atom.get_name()
        self.alt_loc = atom.get_altloc() #Return alternative location specifier.
        self.fullid=atom.get_full_id()
        #('3ukr_test', 0, 'G', (' ', 113, ' '), ('CA', ' '))
        self.res = atom.get_parent().get_resname()
        self.chain = atom.get_full_id()[2]
        self.res_no = int(self.fullid[3][1])
        self.model = atom.get_full_id()[1]
        self.icode = ""
        if atom.is_disordered()==1:
            self.icode = "D"
             # 1 if the residue has disordered atoms
#            self.icode = pdbString[26].strip()#code for insertion residues
#             # Starting co-ordinates of atom.
        self.init_x = atom.get_coord()[0]
        self.init_y = atom.get_coord()[1]
        self.init_z = atom.get_coord()[2]
#             # Current co-ordinates of atom.
        self.x = float(atom.get_coord()[0])
        self.y = float(atom.get_coord()[1])
        self.z = float(atom.get_coord()[2])
#             
        self.occ = atom.get_occupancy()
        self.temp_fac = atom.get_bfactor()
        try:
            self.elem = atom.get_element()
        except:
                self.elem=""
        self.charge=""  
            #Mass of atom as given by atomicMasses global constant. Defaults to 1.
        self.mass = atomicMasses.get(self.atom_name[0])
        if not self.mass:
                self.mass = 1.0
#             # True if atom is the terminal of a chain. Automatically false until modified.
        #vdW of an atom
        try: self.vdw = vdw_radii[self.atom_name[0]]
        except: self.vdw = 1.7
        if not self.vdw: self.vdw = 1.7
        self.isTerm = False
        #for atoms in a grid box, store grid indices
        self.grid_indices = []
    
    def __repr__(self):
        return '('+ self.get_res() +' '+ str(self.res_no) + ' '+self.chain + ': ' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'


    def copy(self):
        """
        
        Return:
            Copy of the Atom instance.
        """
        atom = BioPyAtom([])
        atom.record_name = self.record_name
        atom.serial = self.serial
        atom.atom_name = self.atom_name
        atom.alt_loc = self.alt_loc
        atom.res = self.res
        atom.chain = self.chain
        atom.res_no = self.res_no
        atom.icode = self.icode
        atom.init_x = self.init_x
        atom.init_y = self.init_y
        atom.init_z = self.init_z
        atom.x = self.x
        atom.y = self.y
        atom.z = self.z
        atom.occ =self.occ
        atom.temp_fac = self.temp_fac
        atom.elem = self.elem
        atom.charge = self.charge
        atom.mass = self.mass
        atom.isTerm = self.isTerm
        return atom

    def get_mass(self):
        """
        
        Return:
            Atom mass.
        """
        return self.mass

    def distance_from_init_position(self):
        """
        
        Return:
            Distance from initial position.
        """
        return ((self.x - self.init_x)**2 + (self.y - self.init_y)**2 + (self.z - self.init_z)**2)**0.5

    # Was 'distance_from_atom'
    def distance_from_atom(self, atom):
        """
        
        Return:
            Distance from another atom.
        """
        return ((self.x - atom.x)**2 + (self.y - atom.y)**2 + (self.z - atom.z)**2)**0.5

    def reset_position(self):
        """
        
        Translate atom back in its initial position.
        
        Return:
            initial position co-ordinates of atom.
        """
        self.x = self.init_x
        self.y = self.init_y
        self.z = self.init_z

    def change_init_position(self):
        """
        
        Change initial position co-ordinates of atom to current position.
        
        Return:
            new initial position co-ordinates of atom.
        """
        self.init_x = self.x
        self.init_y = self.y
        self.init_z = self.z

    def translate(self, x, y, z):
        """
        
        Translate the atom.
        
        Arguments:
        
            *x, y, z*
                distance in Angstroms in respective directions to move atom.
        
        Return:
            Translate atom object.
            
        """
        self.x += x
        self.y += y
        self.z += z

    def map_grid_position(self, densMap):
        """
                          
        Arguments:   
            *densMap*
                EM map object consisting the 3D grid of density values.
                
        Return:
              The co-ordinates and density value of the grid point in a density map closest to this atom.
              Return 0 if atom is outside of map.
        """
        x_origin = densMap.x_origin
        y_origin = densMap.y_origin
        z_origin = densMap.z_origin
        apix = densMap.apix
        x_size = densMap.x_size
        y_size = densMap.y_size
        z_size = densMap.z_size
        x_pos = int((self.getX()-x_origin)/apix)
        y_pos = int((self.getY()-y_origin)/apix)
        z_pos = int((self.getZ()-z_origin)/apix)
        if((x_size > x_pos >= 0) and (y_size > y_pos >= 0) and (z_size > z_pos >= 0)):
            return (x_pos, y_pos, z_pos, self.mass)
        else:
            return 0

    def matrix_transform(self, rot_mat):
        """
        
        Transform atom using a 3x3 matrix
                  
        Arguments:   
            *rot_mat*
                a 3x3 matrix instance.
                
        Return:
            Transformed position of atom object
        """
        atom_mat = matrix([[self.x],[self.y],[self.z]])
        new_pos = rot_mat*atom_mat
        self.x = float(new_pos[0])
        self.y = float(new_pos[1])
        self.z = float(new_pos[2])

    def get_pos_vector(self):
        """
        
        Return:
            Vector instance containing 3D coordinates of the atom.
        """
        return Vector.Vector(self.x, self.y, self.z)

    def get_pos_mass(self):
        """
        
        Return:
            An array containing Vector instances containing 3D coordinates of the atom and and its corresponding mass.
        """
        return [self.x, self.y, self.z, self.mass]
        
    def get_x(self):
        """
        
        Return:
            x co-ordinate of atom.
        """
        return float(self.x)
    
    def get_y(self):
        """
        
        Return: 
            y co-ordinate of atom.
        """
        return float(self.y)
    
    def get_z(self):
        """
        
        Return:
            z co-ordinate of atom.
        """
        return float(self.z)
    
    def set_x(self, mod):
        """
        
        Change the x co-ordinate of an atom based on the argument.
        
        Arguments:
            *mod*
                float value
        Return:
            new x co-ordinate
        """
        self.x = mod
    
    def set_y(self, mod):
        """
        
        Change the y co-ordinate of an atom based on the argument.
        
        Arguments:   
            *mod*
                float value
        Return:
            new y co-ordinate
        """
        self.y = mod
    
    def set_z(self, mod):
        """
        
        Change the z co-ordinate of an atom based on the argument.
        
        Arguments:   
            *mod*
                float value
        Return:
            new x co-ordinate
        """
        self.z = mod
    
        
    def get_name(self):
        """
        atom name (ie. 'CA' or 'O')
        
        Return: 
            atom name.
        """
        return self.atom_name
    
    def get_res(self):
        """
        
        Return:
            three letter residue code corresponding to the atom (i.e 'ARG').
        """
        return self.res

    def get_res_no(self):
        """
        
        Return:
            residue number corresponding to the atom.
        """
        return self.res_no

    def get_id_no(self):
        """
        
        Return: 
            string of atom serial number.
        """
        return self.serial

    def _writeTerm(self):
        line = ''
        line += 'TER'.ljust(6)
        line += str(self.serial+1).rjust(5)+' '
        line += ''.center(4)
        line += self.alt_loc.ljust(1)
        line += self.res.ljust(4)
        line += self.chain.ljust(1)
        line += str(self.res_no).rjust(4)
        return line
    

    def write_to_PDB(self):
        """
        
        Writes a PDB ATOM record based in the atom attributes to a file.
        """
        line = ''
        line += self.record_name.ljust(6) 
        line += str(self.serial).rjust(5)+' '
        line += self.atom_name.center(4)
        line += self.alt_loc.ljust(1)
        line += self.res.ljust(3)+' '
        line += self.chain.ljust(1)
        line += str(self.res_no).rjust(4)
        line += str(self.icode).ljust(1)+'   '
        x = '%.3f' % self.x
        y = '%.3f' % self.y
        z = '%.3f' % self.z
        line += x.rjust(8)
        line += y.rjust(8)
        line += z.rjust(8)
        occ = '%.2f'% float(self.occ)
        temp_fac = '%.2f'% float(self.temp_fac)
        line += occ.rjust(6)
        line += temp_fac.rjust(6)+'          '
        line += self.elem.strip().rjust(2)
        line += self.charge.strip().ljust(2)
        return line + '\n'

    #PAP Addition
    def rotate_by_quaternion(self, q_param):
        #print 'inside quaternion rotation for Atom object'
        #create a quaternion object from atom object    
        l = [0.0,self.x,self.y,self.z]
        atom_quat = Quaternion(l)
        q = Quaternion(q_param)
        #create a conjugate of the rotation quaternion
        q_conjugate = q.conjuate(q)
        #perform quaternion rotation
        resultant_quat = q.multiply_3(atom_quat,q_conjugate,q)
        #resultant_quat = q.multiply_3(q,atom_quat,q_conjugate)

        w, x, y, z = resultant_quat.param
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

class _BioPy_Residue:
    """"
    
    A class representing a Residues, use instead the residues option if Biopy
    
    """
    #here to be consistent with the original parser of TEMPY
    #the implemented use of BIOpython is more efficient and retrieve more informations
    """"""
    def __init__(self, res_no, atom_list):
        self.res_no = res_no
        self.atom_list = atom_list[:]

    def _copy(self):
        """
        
        Return:
            Copy of Residues instance.
        """
        newAtomList = []
        for atom in self.atom_list:
            newAtomList.append(atom.copy())
        return _BioPy_Residue(self.res_no, newAtomList)

    def _translate(self, x, y, z):
        """
        
        Translate the Structure.
        
        Arguments:
        
            *x, y, z*
                distance in Angstroms to move structure in respective x, y, and z directions. 
        
        Return:
            Translate Structure object.
            
        """
        for atom in self.atom_list:
            atom.translate(x,y,z)

    
class BioPy_Structure:
    """
    
    A class representing a bjectStructure o, as read from a PDB file using Bio.PDB in Biopython.
    
    
    """

    def __init__(self, atomList, filename='Unknown', header='', footer =''):
        """
        
        Initialise using a string of the relevant pdb file name or a numpy array of Atom objects.
        
        Arguments:
            *pdbFileOrList*
                String of pdb file name or array of Atom objects
            
        """
        self.header = header
        self.footer = footer
        self.filename = filename
        if type(atomList) == ndarray:
            self.atomList = atomList[:]
        if type(atomList) == list:
            self.atomList = array(atomList)
        #Centre of mass calculations
        self.CoM = self.calculate_centre_of_mass()
        ##self.initCoM = self.CoM.copy()
    
    def __getitem__(self, index):
        return self.atomList[index]

    def __len__(self):
        return len(self.atomList)

    def __repr__(self):
        if not self.filename == 'Unknown':
            s =  'Filename: ' + self.filename + '\n'
        else: 
            s = ''
        s += 'No Of Atoms: ' + str(len(self))  + '\n'
        s += 'First Atom: ' + str(self.atomList[0]) + '\n'
        s += 'Last Atom: ' + str(self.atomList[-1]) + '\n'
        return s

    def _write_to_PDB_old(self, filename):
        """
        
        Write Structure instance to PDB file.
        
        Arguments:
        
            *filename*
                output filename.
            
        """
        if filename[-4:] == '.pdb':
            g = open(filename, 'w')
        else:
            g = open(filename+'.pdb', 'w')
        header='''EXPDTA    MODEL GENERATE WITH TEMPY
REMARK    MODEL GENERATE WITH TEMPY
'''
        g.write(header)
        for x in self.atomList:
            g.write(x.write_to_PDB())
            if x.isTerm:
                line = x._writeTerm()+'\n'
                g.write(line)
        g.write(self.footer)
        g.close()

    def write_to_PDB(self, filename,hetatom=False):
        """
        
        Write Structure instance to PDB file.
        
        Arguments:
        
            *filename*
                output filename.
            
        """
        if filename[-4:] == '.pdb':
            g = open(filename, 'w')
        else:
            g = open(filename+'.pdb', 'w')
        header='''EXPDTA    MODEL GENERATE WITH TEMPY
REMARK    MODEL GENERATE WITH TEMPY
'''
        g.write(header)
        structList=self.split_into_chains()
        hetatmlist=[]
        chainlisttot=[]
        last_prot_num=0
        for chain in range(len(structList)):
            chainlist=[]
            for x in structList[chain].atomList:
                if x.record_name=="HETATM":
                    #hetatmlist.append(x.write_to_PDB())
                    hetatmlist.append(x.copy())
                else:
                    chainlist.append(x.copy())
            chainlisttot.append(BioPy_Structure(chainlist))
        for strchain in range(len(chainlisttot)):
            if strchain==0:
                chainlisttot[strchain].renumber_atoms()
            else:
                start_num=chainlisttot[strchain-1].atomList[-1].serial
                #+2 there is ter 
                chainlisttot[strchain].renumber_atoms(start_num=start_num+2)
            for x in chainlisttot[strchain].atomList:
                g.write(x.write_to_PDB())
                last_prot_num+=1
            line = chainlisttot[strchain].atomList[-1]._writeTerm()+'\n'
            g.write(line)
            #last_prot_num+=1

        #print last_prot_num
        if hetatom==True:
            hetstr=BioPy_Structure(hetatmlist)
            hetchain=hetstr.split_into_chains()
            for chain in range(len(hetchain)):
                if chain==0:
                #hetchain[chain].renumber_atoms()
                    #print last_prot_num
                    hetchain[chain].renumber_atoms(start_num=last_prot_num+2)

                else:
                    start_num=hetchain[chain-1].atomList[-1].serial
                #+2 there is ter 
                    hetchain[chain].renumber_atoms(start_num=start_num+2)

                for xhet in hetchain[chain].atomList:
                    g.write(xhet.write_to_PDB())
                line = hetchain[chain].atomList[-1]._writeTerm()+'\n'
                g.write(line)
        lineend = ''
        lineend += 'ENDMDL'.ljust(6)+'\n'
        g.write(lineend)
        g.write(self.footer)
        g.close()

    def copy(self):
        """
        
        Return:
            Copy of Structure instance.
        
        """
        newAtomList = []
        for atom in self.atomList:
            newAtomList.append(atom.copy())
        return BioPy_Structure(newAtomList)

    #===========================================================================
    # def add_het_atoms(self, filename):
    #     """Add het atoms from filename to this Structure instance.
    #     #NB: DO NOT USE IT ATM
    #     """
    #     
    #     self.atomList[-1].isTerm = True
    #     HETATMList = []
    #     f = open(filename, 'r')
    #     for line in f.readlines():
    #         info = line.split()
    #         if(info[0][:6] == 'HETATM'):
    #             HETATMList.append(Atom(line, False))
    #     f.close()
    #     print len(HETATMList)
    #     self.atomList = append(self.atomList, HETATMList)
    #     self.renumber_atoms()
    #===========================================================================

    def calculate_centre_of_mass(self):
        """
        
        Return:    
            Center of mass of structure as a Vector instance.
        
        """
        x_momentTotal = 0.0
        y_momentTotal = 0.0
        z_momentTotal = 0.0
        massTotal = 0.0
        for atom in self.atomList:
            x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            m = atom.get_mass()
            x_momentTotal += x*m
            y_momentTotal += y*m
            z_momentTotal += z*m
            massTotal += m
            x_CoM = x_momentTotal/massTotal
            y_CoM = y_momentTotal/massTotal
            z_CoM = z_momentTotal/massTotal
        return Vector.Vector(x_CoM, y_CoM, z_CoM)

    def translate(self, x, y, z):
        """
        
        Translate the structure.
        
        Arguments:
        
            *x, y, z*
                distance in Angstroms in respective directions to move structure.
                
        Return:
            Translated Structure instance
            
        """
        for atom in self.atomList:
            atom.translate(x,y,z)
        self.CoM = self.CoM.translate(x,y,z)
        

    def rotate_by_axis_angle(self, x, y, z, turn, x_trans=0, y_trans=0, z_trans=0, rad=False, com=False):
        """
        
        Rotate the Structure instance around its centre.
        
        Arguments:
        
            *turn*
                angle (in radians if rad == True, else in degrees) to rotate map.
            *x,y,z*
                axis to rotate about, ie. x,y,z =  0,0,1 rotates the structure round the xy-plane.
            *x_trans, y_trans, z_trans*
                extra translational movement if required.
            *com*
                centre of mass around which to rotate the structure. If False, rotates around centre of mass of structure.
                
        """

        mat = Vector.axis_angle_to_matrix(x, y, z, turn, rad)

        if not com:
            com = self.CoM.copy()

        newcom = com.matrix_transform(mat)
        offset = com-newcom
        self.matrix_transform(mat)
        self.translate(x_trans+offset.x, y_trans+offset.y, z_trans+offset.z)
        

    def rotate_by_euler(self, x_turn, y_turn, z_turn, x_trans=0, y_trans=0, z_trans=0, rad=False, com=False):
        """
        
        Rotate this Structure instance around its centre.
        
        Arguments:
            *x_turn,y_turn,z_turn*
                Euler angles (in radians if rad == True, else in degrees) used to rotate structure, in order XYZ.
            *x_trans, y_trans, z_trans*
                extra translational movement if required.
            *com*
                centre of mass around which to rotate the structure. If False, rotates around centre of mass of structure.
                
        """

        mat = Vector.euler_to_matrix(x_turn, y_turn, z_turn, rad)

        if not com:
            com = self.CoM.copy()

        newcom = com.matrix_transform(mat)
        offset = com-newcom
        self.matrix_transform(mat)
        self.translate(x_trans+offset.x, y_trans+offset.y, z_trans+offset.z)

    #PAP addition
    def rotate_by_quaternion(self, q, com=False):
        #print 'inside quaternion rotation for Structure object'
        if not com:
            com = self.CoM.copy()
        self.translate(-com.x, -com.y, -com.z)
        for atom in self.atomList:
                atom.rotate_by_quaternion(q)
        self.translate(com.x, com.y, com.z)

    def randomise_position(self, max_trans, max_rot, v_grain=30, rad=False, verbose=False):
        """
        
        Randomise the position of the Structure instance using random rotations and translations. 
                  
        Arguments:   
            *max_trans*
                Maximum translation permitted
            *max_rot*
                Maximum rotation permitted (in degree if rad=False)
            *v_grain*
                Graning Level for the generation of random vetors (default=30)
        Return:
            Transformed position of Structure object
        """
        t_v = Vector.random_vector(-v_grain, v_grain).unit()
        r_v = Vector.random_vector(-v_grain, v_grain).unit()

        if max_trans <= 0:
            t_dist = 0
        else:
            t_dist = randrange(max_trans)

        if max_rot <= 0:
            r_ang = 0
        else:
            r_ang = randrange(max_rot)
        t_v = t_v.times(t_dist)
        self.rotate_by_axis_angle(r_v.x, r_v.y, r_v.z, r_ang, t_v.x, t_v.y, t_v.z, rad=rad)
        #self.translate(t_v.x, t_v.y, t_v.z)
        if verbose:
            print (r_v.x, r_v.y, r_v.z, r_ang, t_v.x, t_v.y, t_v.z)
        return (r_v.x, r_v.y, r_v.z, r_ang, t_v.x, t_v.y, t_v.z)

    def matrix_transform(self, matrix):
        """
        
        Transform Structure using a 3x3 transformation matrix
                  
        Arguments:   
            *rot_mat*
                a 3x3 matrix instance.
                
        Return:
            Transformed position of Structure object
        """
        for atom in self.atomList:
            atom.matrix_transform(matrix)
        self.CoM = self.CoM.matrix_transform(matrix)

    def reorder_residues(self):
        """
        
        Order residues in atom list by residue number. 
        (NOTE: Does not check for chain information - split by chain first).
        
        """
        self.atomList = list(self.atomList)
        self.atomList.sort(cmp=lambda x,y: int(x.res_no)-int(y.res_no))
        self.atomList = array(self.atomList)
        
    def renumber_residues(self, startRes=1, missingRes=[]):
        """
        
        Renumber the structure starting from startRes.
        Missing number list to add.
        
        Arguments:
            *startRes*
                Starting residue number for renumbering
            *missingRes*
                A list of missing residue numbers to add 
                 
        """
        resNo = startRes
        currentRes = self.atomList[0].res_no
        for x in self.atomList:
            if x.res_no == currentRes:
                x.res_no = resNo
            else:
                currentRes = x.res_no
                resNo +=1
                while resNo in missingRes:
                    resNo +=1
                x.res_no = resNo

    def renumber_atoms(self,start_num=1):
        """
        
        Renumber the atoms in the structure. 
        After renumbering the starting atom number will be 1 unless start_num
        
        """
        for x in range(len(self.atomList)):
            if (x+1) < 99999:
                self.atomList[x].serial = x+start_num
            else:
                self.atomList[x].serial = '*****'

    def rename_chains(self, chain_list=False):
        """
        
        Rename chain name based on the list of new chain names
        
        Arguments:
             *chain_list*
                 List of chain names
                 If False rename in alphabetical order.

        """
        
        if not chain_list:
            chain_list = alphabet
        else:
            noc = self.no_of_chains()
            if len(chain_list) != self.no_of_chains():
                print('No. of chains in structure = '+str(noc))
                print('Length of chain list = '+str(len(chain_list)))
                print('Chains not changed.')
                return
        ch = self.atomList[0].chain
        renum = 0
        for atom in self.atomList:
            if atom.chain == ch:
                atom.chain = chain_list[renum]
            else:
                renum += 1
                ch = atom.chain[:]
                atom.chain = chain_list[renum]

    def split_into_chains(self):
        """
         
         Split the structure into separate chains and returns the list of Structure instance for each chain. 
         
        """
        structList = []
        currentChain = self.atomList[0].chain
        currentStruct = []
        for x in self.atomList:
            if x.chain == currentChain:
                currentStruct.append(x.copy())
            else:
                currentChain = x.chain
                structList.append(BioPy_Structure(currentStruct))
                currentStruct = [x.copy()]
        structList.append(BioPy_Structure(currentStruct))
        return structList

    def no_of_chains(self):
        """
        Return:
            the number of chains in the Structure object
        """
        a = self.split_into_chains()
        return len(a)
                
    def reset_position(self):
        """
        
        Translate structure back into initial position.
        
        """
        for x in self.atomList:
            x.reset_position()
        self.CoM = self.initCoM.copy()

    def change_init_position(self):
        """
        
        Change initial position of structure to current position.
        
        """
        for atom in self.atomList:
            atom.change_init_position()
        self.init_CoM = self.CoM.copy()
    
    def RMSD_from_init_position(self, CA=False):
        """
        
        Return RMSD of structure from initial position after translation.
        
        Arguments:
            *CA*
                True will consider only CA atoms.
                False will consider all atoms.
        Return:
            RMSD in angstrom
            
        """
        dists = []
        for x in self.atomList:
            if CA:
                if x.atom_name == 'CA':
                    dists.append(x.distance_from_init_position())
            else:
                dists.append(x.distance_from_init_position())
        dists = array(dists)
        return dists.mean()

    def RMSD_from_same_structure(self, otherStruct, CA=False, write=True):
        """
        
        Return the RMSD between two structure instances.
        
        Arguments:
            *otherStruct*
                Structure instance to compare, containing the same number of atoms as the target instance.
            *CA*
                True will consider only CA atoms.
                False will consider all atoms.
        Return:
            RMSD in angstrom
 
        """            
        dists = []
        for a in range(len(self.atomList)):
            if CA:
                if self.atomList[a].atom_name == 'CA':

                        if otherStruct.atomList[a].atom_name == 'CA':
                            dists.append(self.atomList[a].distance_from_atom(otherStruct.atomList[a]))
                            print(self.atomList[a].atom_name, otherStruct.atomList[a].atom_name,self.atomList[a].res_no, otherStruct.atomList[a].res_no)
 
                        else:
                            pass
            else:
                dists.append(self.atomList[a].distance_from_atom(otherStruct.atomList[a]))
        dists = array(dists)
        return dists.mean()

    
    def get_vector_list(self):
        """
        
        Return:
            Array containing 3D Vector instances of positions of all atoms.
        
        """
        v = []
        for atom in self.atomList:
            v.append(atom.get_pos_vector())
        return array(v)

    def get_pos_mass_list(self):
        """
        
        Return:
            Array containing Vector instances of positions of all atoms and mass.
        
        """
        v = []
        for atom in self.atomList:
            v.append(atom.get_pos_mass())
        return array(v)

    def get_extreme_values(self):
        """
        
        Return:
            A 6-tuple containing the minimum and maximum of x, y and z co-ordinates of the structure.
            Given in order (min_x, max_x, min_y, max_y, min_z, max_z).
        
        """
        min_x = self.atomList[0].get_x()
        max_x = self.atomList[0].get_x()
        min_y = self.atomList[0].get_y()
        max_y = self.atomList[0].get_y()
        min_z = self.atomList[0].get_z()
        max_z = self.atomList[0].get_z()
        for atom in self.atomList[1:]:
            if atom.get_x() < min_x:
                min_x = atom.get_x()
            if atom.get_x() > max_x:
                max_x = atom.get_x()
            if atom.get_y() < min_y:
                min_y = atom.get_y()
            if atom.get_y() > max_y:
                max_y = atom.get_y()
            if atom.get_z() < min_z:
                min_z = atom.get_z()
            if atom.get_z() > max_z:
                max_z = atom.get_z()
        return (min_x, max_x, min_y, max_y, min_z, max_z)

                
    def get_atom_list(self):
        """
        
        Return:
            An array containing Atom instances of positions of all atoms as:
            [(RES 1 A: x,y,z), ... ,(RES2 1 A: x1,y1,z1)].
        
        """
        alist = []
        for x in self.atomList:
            alist.append(x.copy())
        return alist

    def find_same_atom(self, atom_index, otherStruct):
        """
        Find if an atom exists in the compared structure, based on atom index.
 
        Arguments:
            *atom_index*
                atom number
            *otherStruct*
                a structure object to compare
  
        Return:
            If a match is found, it returns the atom object; else it returns a string reporting the mismatch.
        
        """
        atom = self.atomList[atom_index]
        for x in otherStruct.atomList:
            if x.res_no == atom.res_no and x.atom_name == atom.atom_name and atom.res == atom.res and x.chain == atom.chain:
                #print x
                return x
        return "No match of atom index %s in structure %s"%(atom_index,otherStruct)

    def get_chain_list(self):
        """
        
        Return:
            A list of chain ID.
        
        """
        chain_list=[]
        for x in self.atomList:
            if x.chain not in chain_list:
                chain_list.append(x.chain) 
        return chain_list
    
    def get_chain(self, chainID):
        """
        
        Return:
            New Structure instance with only the requested chain.
        
        """
        newAtomList = []
        for x in self.atomList:
            if x.chain == chainID:
                newAtomList.append(x.copy())
        if len(newAtomList)!=0:
            return BioPy_Structure(newAtomList)
        else:
            print("Warning no chain " + chainID + "found")
    
    
    def get_selection(self, startRes, finishRes, chain=''):
        """
        
        Get a new Structure instance for the selected residues range without considering residues chain. 
        
        Arguments:
            *startRes*
                Start residue number
            *finishRes*
                End residue number 
                
        Return:
            New Structure instance
        
        """
        newAtomList = []
        for x in self.atomList:
            cur_chain = x.chain
            if chain:
                if not cur_chain == chain: continue
            if(x.get_res_no() >= int(startRes) and x.get_res_no() <= int(finishRes)):
                newAtomList.append(x.copy())
        return BioPy_Structure(newAtomList)
    
    

    def break_into_segments(self, rigid_list,chain=''): 
        """
        
        Return a list of Structure instance based on the rigid body list.
        
        Arguments:
        
            *rigid list*
                list of rigid body defined as:
            
                    [[riA,rfA],..,[riB,rfB]]
            
                where :
            
                    riA is the starting residues number of segment A.
                    rfA is the final residues number of segment A.
        Return:
            List of TEMPy Structure instance

            
        """
        structList = []
        for r in rigid_list:
            fstStruct = self.get_selection(r[0], r[1],chain)
            nxtStructs = []
            for x in range(2, len(r), 2):
                nxtStructs.append(self.get_selection(r[x], r[x+1],chain))
            #print nxtStructs
            if len(nxtStructs) != 0:
                structList.append(fstStruct.combine_structures(nxtStructs))
            else:
                structList.append(fstStruct)
        if len(structList) ==0:
        	print("Error: Residues not in PDB.")
        	sys.exit()
        else:
        	return structList
    
    #PAP ADDITION       
    def get_chain_ca(self, chainID):
        newAtomList = []
        for x in self.atomList:
            if x.chain == chainID and x.atom_name == "CA":
                newAtomList.append(x.copy())
        if len(newAtomList)!=0:
            return BioPy_Structure(newAtomList)
        else:
            print("Warning no chain " + chainID + " found")

    #PAP ADDITION
    def get_rgyration(self):
        vc = self.calculate_centre_of_mass()
        rg = 0.0
        num = 0.0
        den = 0.0
        for x in self.atomList:
                vx = x.get_pos_vector()
                dist = vx.dist(vc)
                num = num + x.get_mass() * (dist*dist)
                den = den + x.get_mass()
        rg = sqrt(num/den)
        return rg
    #PAP ADDITION
    def calculate_centroid(self):
        """Return centre of mass of structure as a Vector instance."""
        x_Total = 0.0
        y_Total = 0.0
        z_Total = 0.0
        natom = 0.0
        for atom in self.atomList:
            x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            x_Total += x
            y_Total += y
            z_Total += z
            natom += 1
            x_CoM = x_Total/natom
            y_CoM = y_Total/natom
            z_CoM = z_Total/natom
        return Vector.Vector(x_CoM, y_CoM, z_CoM)
    
#===============================================================================
# #add IF 19/2/2013 of residues number
# #we decided to not have this.
#     def _get_RB(self, rigid_list):
#         """
#         
#         Return a Structure instance of the selected rigid body of non consecutively segment from a Structure instance:
#         
#         Arguments:
#             *rigid list* 
#                 list of rigid body defined as:
#             
#                 [riA,rfA,..,riB,rfB]
#             
#                 where :
#             
#                 riA= starting residues of segment A
#                 rfA=final residues of segment A
#             
#                 so that Rigid body Structure is formed by non consecutively segment.
#                 
#         """
#         structList = []
#         nxtStructs=[]
#         fstStruct = self.get_selection(rigid_list[0],rigid_list[1])         
#         if len(rigid_list)>2:
#             for r in range(2, len(rigid_list), 2):
#                 nxtStructs.append(self.get_selection(rigid_list[r], rigid_list[r+1]))
#                 if len(nxtStructs) != 0:
#                     structList.append(fstStruct.combine_structures(nxtStructs))
#                 else:
#                     structList.append(fstStruct)
#     #return the last structureList that contains only the atoms of the defined RB
#         return structList[-1]
#===============================================================================

    def combine_structures(self, structList):
        """
        
        Add a list of Structure instance to the existing structure.
        
        Arguments:
            *structList*
                list of Structure instance
        Return:
            New Structure Instance
        
        """
        atomList = self.atomList.copy()
        for s in structList:
            atomList = append(atomList, s.atomList)
        return BioPy_Structure(atomList)
    
    def combine_SSE_structures(self, structList):
        """
        
        Combine a list of Structure instance into one and return a new Structure instance.
        
        Arguments:
            *structList*
                list of Structure instance
        Return:
            New Structure Instance
        
        """
        atomList =[]
        for s in structList:
            atomList = append(atomList, s.atomList)
        return BioPy_Structure(atomList)
    
    def get_selection_more_than(self, startRes):
        """
        
        Get a Structure instance comprising all residues with their residue numbers greater than startRes.
        
        Arguments:
            *startRes*
                a residue number
        
        Return:
            A Structure instance
        """
        newAtomList = []
        for x in self.atomList:
            if(x.get_res_no() >= int(startRes)):
                newAtomList.append(x.copy())
        return BioPy_Structure(newAtomList)
            
    def get_selection_less_than(self, endRes):
        """
        
        Get a Structure instance comprising all residues with their residue numbers less than endRes.
        
        Arguments:
            *endRes*
                a residue number
        
        Return:
            A Structure instance
        
        """
        newAtomList = []
        for x in self.atomList:
            if(x.get_res_no() <= endRes):
                newAtomList.append(x.copy())
        return BioPy_Structure(newAtomList)

    def get_residue(self, resNo):
        """
       
       Get the residue corresponding to the residue number.
        
        Arguments:
            *resNo*
                Residues number
        Return: 
            Returns a Residues instance. 

        """
        return BioPy_Structure([x.copy() for x in self.atomList if x.get_res_no() == int(resNo)])
    
    def get_atom(self, index):
        """
        
        Return specific atom in Structure instance.
        
        Arguments:
            *index*
                Index of the atom
        Return: 
            Returns an Atom instance. 

        """
        return self.atomList[int(index)]

    def get_backbone(self):
        """
        
        Return:
            Structure instance with only the backbone atoms in structure.
        
        """
        backboneList = []
        #print self.atomList
        for atom in self.atomList:#error again due to aold atom list 
            if(atom.get_name() == 'CA' or atom.get_name() == 'N' or atom.get_name() == 'C'):
                backboneList.append(atom.copy())
        return BioPy_Structure(backboneList[:])

    def get_CAonly(self):
        """
        
        Return:
            Structure instance with only the backbone atoms in structure.
        
        """
        backboneList = []
        #print self.atomList
        for atom in self.atomList:#error again due to aold atom list 
            if(atom.get_name() == 'CA'):
                backboneList.append(atom.copy())
            else:
                pass
        return BioPy_Structure(backboneList)

    def vectorise(self):
        vectorList = []
        vList = []
        for x in self.atomList:
            vectorList.append(x.get_pos_vector())
        for y in range(len(vectorList)-1):
            vList.append(vectorList[y]-(vectorList[y+1]))
        return vList
        
    def get_torsion_angles(self):
        """
        
        Return:
            List of torsion angles in Structure instance.
        """
        vectorList = self.vectorise()
        angles = []
        for v in range(len(vectorList)-2):
            angles.append(Vector.altTorsion(vectorList[v], vectorList[v+1].reverse(), vectorList[v+2]))
        return angles


    def get_prot_mass_from_res(self,Termini=False):
        # ADD by IF 22-4-2013 
        #from Harpal  code calculation of mass from seq.
        """
        Calculates Mass (kDa) of the Structure instance, from average mass
        
        """      
        #NOTE: problem with Residues class from Structure class need a more elegand way of doing it
        #PROBLEM ONLY 1 CHAIN READ
        #atm seq_list_resno to create the sequence from pdb
            
        aa = {'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E','SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','SEC':'U','GLY':'G','PRO':'P','ALA':'A','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V'}
        
        #based on http://web.expasy.org/findmod/findmod_masses.html
        mass_tot=0
        str=self.copy()
        seq_string=''
        for chain in str.split_into_chains():
            seq_list_resno=[]
            seq_str_aa=[]
            for x in chain.atomList:
                if x.res in aa.keys():
                    if x.res_no not in seq_list_resno: 
                        seq_list_resno.append(x.res_no)
                        if x.res not in aa.keys():
                            seq_string+="x"
                        res_singleletter=aa[x.res]
                        seq_string+="%s"%res_singleletter
                        mass_tot += aa_mass[res_singleletter]
                else:
                    pass
        if Termini:
            mass_tot += 17.992 
 
        return float(mass_tot/1000)
     
    def get_prot_mass_from_atoms(self):
        """
        Calculates Mass (kDa) of the Structure instance, from average mass. Atoms based
        use get_prot_mass_from_res is more accurate.
        """      

        # ADD by IF 22-4-2013 
        #problem with this are the PDBs not cleaned with ANISU or similar problems
        mass_tot=0
        for x in self.atomList:
            mass_tot+=x.get_mass()
        return float(mass_tot/1000)
    
    
    
#PRIVATE 
#functions in development for symmetry operations.
    
#===============================================================================
#     
#     def build_C_symmetric_model_around_line(self, l, c, noOfUnits):
#         output = [self.copy()]
#         output[0].rename_chains()
#         angle = 360./noOfUnits
#         
#         for x in range(1,noOfUnits):
#             self.rotate_by_axis_angle(l.x, l.y, l.z, angle*x, com=c)
#             next_unit = self.copy()
#             next_unit.change_init_position()
#             for atom in next_unit.atomList:
#                 atom.chain = alphabet[x]
#             output.append(next_unit)
#             self.reset_position()
# 
#         output = output[0].combine_structures(output[1:])
#         output.CoM = c.copy()
#         output.renumber_atoms()
#         return output
# 
#     def build_C_symmetric_model(self, symm_axis, noOfUnits, x_point, y_point, z_point):
#         c = Vector.Vector(x_point, y_point, z_point)
#         if symm_axis == 'x':
#             l = Vector.Vector(1,0,0)
#         elif symm_axis == 'y':
#             l = Vector.Vector(0,1,0)
#         else:
#             l = Vector.Vector(0,0,1)
#         return self.build_C_symmetric_model_around_line(l, c, noOfUnits)
# 
#     def build_C_symmetric_model_using_map(self, noOfUnits, symm_axis, densMap):
#         centre = densMap.centre()
#         return self.build_C_symmetric_model(symm_axis, noOfUnits, centre[0], centre[1], centre[2])
#===============================================================================
