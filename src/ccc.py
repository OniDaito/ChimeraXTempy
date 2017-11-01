# vim: set expandtab shiftwidth=2 softtabstop=2:
# ccc
#   "score" - compute the ccc score using tempy

from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom

from .util import chimera_to_tempy_atom, chimera_to_tempy_map

import math

def score(session, atomic_model, map_model, rez):
  ''' Perform the CCC score. Takes a session, a single model and map.'''

  print("Calculating CCC Score")
 
  # make class instances for density simulation (blurring), scoring and plot scores
  blurrer = StructureBlurrer()
  scorer = ScoringFunctions()

  atomlist = []
  for atom in atomic_model.atoms:
    atomlist.append(chimera_to_tempy_atom(atom, len(atomlist)))
    
  bio_atom_structure = BioPy_Structure(atomlist)
  bio_map_structure = chimera_to_tempy_map(map_model)
  map_probe = blurrer.gaussian_blur(bio_atom_structure, rez, densMap=bio_map_structure)
  score = scorer.CCC(bio_map_structure, map_probe)
  print(score)
  return score
