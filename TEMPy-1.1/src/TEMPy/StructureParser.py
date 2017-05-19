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

from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom
import urllib
from numpy import append
import sys


class mmCIFParser:
    """A class to read mmCIF files either directly from the mmCIF or a structure instance from Biopython"""
    def __init__(self):
        pass

    @staticmethod
    def read_mmCIF_file(structure_id, filename,hetatm=False,water= False):
        """
        
        Read mmCIF file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of mmCIF file       
            *filename*
                name of mmCIF file
            *hetatm*
                Boolean representing whether the mmCIF file contains hetatom.
                Default and recommended is False.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.
        
        Return:
            Structure Instance
        """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        p=MMCIFParserBiopy()#permissive default True
        structure=p.get_structure(structure_id, filename)
        return mmCIFParser._biommCIF_strcuture_to_TEMpy(filename,structure,hetatm,water)

    @staticmethod
    def fetch_mmCIF(structure_id, filename,hetatm=False,water= False):
        
        """
        
        Fetch mmCIF file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of mmCIF file       
            *filename*
                name of mmCIF file
            *hetatm*
                Boolean representing whether the mmCIF file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.
        
        Return:
            Structure Instance
         """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        
        p=MMCIFParserBiopy()
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % structure_id
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename)
        return mmCIFParser._biommCIF_strcuture_to_TEMpy(filename,structure,hetatm,water)

    @staticmethod
    def _biommCIF_strcuture_to_TEMpy(filename,structure,hetatm=False,water= False):
            #imported if and when the function is executed.
        """
        PRIVATE FUNCTION to convert to Structure Instance
        filename = name of mmCIF file
        hetatm = Boolean representing whether to add hetatm to the structure.Default and Raccomanded is False.
        water = Boolean representing whether to add water to the structure.Default and Raccomanded is False.
        """
        from Bio.PDB import MMCIFParser as MMCIFParserBiopy
        
        p=MMCIFParserBiopy()
        
        atomList = []
        hetatomList=[]
        wateratomList=[]
        footer = ''
        header = ''
        cif_code=filename.split("/")[-1]#use os.1FAT.cif
        structure_id="%s" % cif_code[:-4]
        structure=p.get_structure(structure_id, filename)
        residues = structure.get_residues()
        for res in residues:
            hetfield=res.get_id()[0]
            if hetfield[0]=="H":
                for atom in res:
                    BioPyAtom(atom)
                    hetatomList.append(BioPyAtom(atom))
            elif hetfield[0]=="W":
                for atom in res:
                    BioPyAtom(atom)
                    wateratomList.append(BioPyAtom(atom))
            else:
                for atom in res:
                    BioPyAtom(atom)
                    atomList.append(BioPyAtom(atom))
        if hetatm:
            atomList = append(atomList, hetatomList)
        if water:
            atomList = append(atomList, wateratomList)
        
        return BioPy_Structure(atomList, filename=filename, header=header, footer=footer)


class PDBParser:
    """A class to read PDB files either directly from the pdb or a structure instance from Biopython"""
    def __init__(self):
        pass

    @staticmethod
    def read_PDB_file(structure_id, filename,hetatm=False,water= False,chain=None):
        """
        
        Read PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
        
        p=PDBParserBiopy(QUIET=True)#permissive default True
#        try:
        structure=p.get_structure(structure_id, filename)
        #except AssertionError:
         #   sys.stderr.write('unknown element in BioPyton\n')
          #  sys.stderr.write(structure_id)
           # sys.exit()
        return PDBParser._bio_strcuture_to_TEMpy(filename,structure,hetatm,water)


    @staticmethod
    def read_PDB_file_BioPy(structure_id, filename,hetatm=False,water= False,chain=None):
        """
        
        Read PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
        
        p=PDBParserBiopy(QUIET=True)#permissive default True
#        try:
        structure=p.get_structure(structure_id, filename)
        #except AssertionError:
         #   sys.stderr.write('unknown element in BioPyton\n')
          #  sys.stderr.write(structure_id)
           # sys.exit()
        return structure
    
    @staticmethod
    def fetch_PDB(structure_id, filename,hetatm=False,water= False):       
        """
 
        Fetch PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
 
        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % structure_id
        p=PDBParserBiopy(QUIET=True)#permissive default True
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename)
        return PDBParser._bio_strcuture_to_TEMpy(filename,structure,hetatm,water)
        
    @staticmethod
    def fetch_PDB_BioPy(structure_id, filename,hetatm=False,water= False):       
        """
 
        Fetch PDB file and create Structure instance based upon it.
           
        Argument:
            *structure_id*
                structure_id code of pdb file       
            *filename*
                name of pdb file
            *hetatm*
                Boolean representing whether the PDB file contains hetatom.
            *water*
               Boolean representing whether to add water to the structure.
               Default and recommended is False.

        Return:
            Structure Instance
        """
        from Bio.PDB import PDBParser as PDBParserBiopy
 
        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % structure_id
        p=PDBParserBiopy(QUIET=True)#permissive default True
        urllib.urlretrieve(url, filename)
        structure=p.get_structure(structure_id, filename,hetatm=hetatm,water=water)
        return structure

    @staticmethod
    def _bio_strcuture_to_TEMpy(filename,structure,hetatm=False,water= False):
            #imported if and when the function is executed.
        """
        PRIVATE FUNCTION to convert to Structure Instance
        filename = name of mmCIF file
        hetatm = Boolean representing whether to add hetatm to the structure.Default and Raccomanded is False.
        water = Boolean representing whether to add water to the structure.Default and Raccomanded is False.
        """
        #from Bio.PDB import PDBParser as PDBParserBiopy
        atomList = []
        hetatomList=[]
        wateratomList=[]
        footer = ''
        header = ''
        #pdb_code=filename.split("/")[-1]#use os.
        #p=PDBParserBiopy()#permissive default True
        #structure_id="%s" % pdb_code[:-4]
        #structure=p.get_structure(structure_id, filename)
        residues = structure.get_residues()
        for res in residues:
            hetfield=res.get_id()[0]
            if hetfield[0]=="H":
                for atom in res:
                    BioPyAtom(atom)
                    hetatomList.append(BioPyAtom(atom))
            elif hetfield[0]=="W":
                for atom in res:
                    BioPyAtom(atom)
                    wateratomList.append(BioPyAtom(atom))
            else:
                for atom in res:
                    BioPyAtom(atom)
                    atomList.append(BioPyAtom(atom))
        if hetatm:
            atomList = append(atomList, hetatomList)
        if water:
            atomList = append(atomList, wateratomList)
        
        return BioPy_Structure(atomList, filename=filename, header=header, footer=footer)
    
    @staticmethod
    def calc_SA(self,pdbfile,rsa=True,outsafile=None):
        assert os.path.isfile(pdbfile)
        if outsafile is None: outsafile = os.path.basename(pdbfile)+'_sa.out'
        #o = open(outsafile,'w')
        cmd = "~/data/packages/freesasa/freesasa-1.1/src/freesasa %s --rsa_file=%s\
         --no-log --radii=naccess"%(pdbfile,outsafile)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        
        
        

