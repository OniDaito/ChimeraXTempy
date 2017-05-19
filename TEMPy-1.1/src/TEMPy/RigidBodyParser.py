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

class RBParser:
	"""A class to read Rigid Body files"""
	def __init__(self):
		pass

	@staticmethod
	def read_FlexEM_RIBFIND_files(file_in,structure_instance,list_out=True):
		"""
		Read a rigid body files in Flex-EM format (text file) using residue numbers in a structure instance.
		Each line describes one rigid body by specifying the initial and final residue of each of the segments in that rigid body 
		(eg, '2 6 28 30' means that residues 2-6 and 28-30 will be included in the same rigid body). 
		We recommend to use the RIBFIND server for accurately identifying Rigid Bodies in a protein structures.
		
		Arguments:
			*file_in*
				Rigid Body File in Flex EM format
			*structure_instance*
				Structure Instance to manipulate
			*list_out*
				True return a list of the Rigid Bodies structure instances (each line in the file).
				False will print them separately.
				
		
		"""
		ssefile = open(file_in, 'rU')
		nsse = 0
		RB_structureinstance_tot=[]
		for line in ssefile:
			if line.startswith("#"):
				pass
			else:
				sselist = []
				nsse += 1
				if len(line) < 3: continue
				tokens = line.split()
				if len(tokens) < 2: continue
				chainID = ''
				for i in range(len(tokens)/2):
					start = int(tokens[i*2].split(':')[0])
					end = int(tokens[i*2+1].split(':')[0])
					sselist.append([start,end])
					if ':' in tokens[i*2]:
						chainID = tokens[i*2].split(':')[1]
						if not tokens[i*2+1].split(':')[1] == chainID:
							print 'Check chain IDs in rigid body file', tokens[i*2]
		#Combine SSEs into one structure
				sse_struct_list = structure_instance.break_into_segments(sselist,chainID)
				sse_struct = structure_instance.combine_SSE_structures(sse_struct_list)
				RB_structureinstance_tot.append(sse_struct.copy())
		
		ssefile.close()
		if list_out:
			try:
				return RB_structureinstance_tot
			except UnboundLocalError:
				print "wrong residues number"
	
		else:
			for structure in RB_structureinstance_tot:
				try:
					return structure
				except UnboundLocalError:
					print "wrong residues number"


	@staticmethod
	def RBfileToRBlist(file_in):
		"""
		Read a rigid body files in Flex-EM format (text file) using residue numbers in a list of segments.
		Each line describes one rigid body by specifying the initial and final residue of each of the segments in that rigid body 
		(eg, '2 6 28 30' means that residues 2-6 and 28-30 will be included in the same rigid body). 
		We recommend to use the RIBFIND server for accurately identifying Rigid Bodies in a protein structures.
		
		Arguments:
			*file_in*
				Rigid Body File in Flex EM format
			*list_out*
				Return a list of segments (each line in the file).
				The list of rigid body is defined as:
			
					[[riA,rfA],..,[riB,rfB]]
			
				where :
			
					riA is the starting residues number of segment A.
					rfA is the final residues number of segment A.
				
		
		"""

		ssefile = open(file_in, 'rU')
		nsse = 0
		RB_structureinstance_tot=[]
		for line in ssefile:
			if line.startswith("#"):
				pass
			else:
				sselist = []
				nsse += 1
				if len(line) < 3: continue
				tokens = line.split()
				if len(tokens) < 2: continue
				for i in range(len(tokens)/2):
					start = tokens[i*2]
					end = tokens[i*2+1]
					sselist.append([start,end])

				RB_structureinstance_tot.append(sselist)
		
		ssefile.close()

		return RB_structureinstance_tot
