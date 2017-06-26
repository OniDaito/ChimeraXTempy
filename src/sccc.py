# vim: set expandtab shiftwidth=2 softtabstop=2:

# sccc
#   "score" - compute the sccc score using tempy
#

from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom

from .util import chimera_to_tempy_atom, chimera_to_tempy_map

import math

def score(session, atomic_model, map_model, rigid_filename, colour_atoms=True):
    print("Calculating SCCC Score")
    rez = 10.0

   
    # the sigma factor determines the width of the Gaussian distribution used to describe each atom
    # TODO - this needs to be a slider or textbox
    sim_sigma_coeff = 0.187
 
    # TODO - replace this with the current loaded map
    #emmap = MapParser.readMRC(map_model_filename)

    # make class instances for density simulation (blurring), scoring and plot scores
    blurrer = StructureBlurrer()
    scorer = ScoringFunctions()

    # read map file
    #emmap=MapParser.readMRC(m)
    
    # read PDB file
    # TODO - we need to change this to current selected PDB Model
    #structure_instance=PDBParser.read_PDB_file('pdbfile',p,hetatm=False,water=False)
  
    atomlist = []

    for atom in atomic_model.atoms:
      atomlist.append(chimera_to_tempy_atom(atom, len(atomlist)))

    bio_atom_structure = BioPy_Structure(atomlist)
    bio_map_structure = chimera_to_tempy_map(map_model)
    SCCC_list_structure_instance=[]
    # read rigid body file and generate structure instances for each segment
    listRB = RBParser.read_FlexEM_RIBFIND_files(rigid_filename, bio_atom_structure)
    
    # score each rigid body segment
    listsc_sccc = []
    print('calculating SCCC')

    for RB in listRB:
      # sccc score

      score_SCCC=scorer.SCCC(bio_map_structure, rez, sim_sigma_coeff, bio_atom_structure, RB, c_mode=False)
      SCCC_list_structure_instance.append(score_SCCC)

      print ('>>', score_SCCC)
      listsc_sccc.append(score_SCCC)
      
      # Colour the atoms based on the rating from white (1.0) to red (0.0)
      # TODO - maybe a faster way? Also 'all_atoms' mentioned in the API doesnt exist but atoms does! :S
      if colour_atoms:
        dr = 255
        dg = 255
        db = 255
        if score_SCCC >= 0.5:
          dr = 255 - int(math.floor(255 * ((score_SCCC - 0.5) * 2.0) ))
          dg = dr
        else:
          db = int(math.floor(255 * (score_SCCC * 2.0)))
          dg = db
 

        residues = []
        for a in RB.atomList:
          if a.res_no not in residues:
            residues.append(a.res_no)

        for r in residues:
          cr = atomic_model.residues[r]
          for catm in cr.atoms:
            catm.color = [dr,dg,db,255] 
          cr.ribbon_color = [dr,dg,db,255]
