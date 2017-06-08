# vim: set expandtab shiftwidth=2 softtabstop=2:

# test_cmd 
#   based on the python3 unittest boilerplate code
#   https://docs.python.org/3.5/library/test.html
#

import unittest

# This might be a bit evil but we need it for now
import sys, os
sys.path.insert(0, "../TEMPy")
 
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ShowPlot import Plot
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser

class TestSCCC(unittest.TestCase):
  # Only use setUp() and tearDown() if necessary
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_tempy_sccc(self):
    ''' Test the tempy sccc score based on the files
    provided. Use this as a baseline for the second
    chimeraX test. '''

    # the sigma factor determines the width of the Gaussian distribution used to describe each atom
    sim_sigma_coeff = 0.187

    path_test = "./"
    m = os.path.join(path_test,'1akeA_10A.mrc')
    p = os.path.join(path_test,'final1_mdcg.pdb')
    r = 10.0
    rb_file = os.path.join(path_test,'rigid_RF.txt')
    
    # make class instances for density simulation (blurring), scoring and plot scores
    blurrer = StructureBlurrer()
    scorer = ScoringFunctions()

    # read map file
    emmap=MapParser.readMRC(m)
    # read PDB file
    structure_instance=PDBParser.read_PDB_file('pdbfile', p, hetatm=False, water=False)
    SCCC_list_structure_instance=[]
    
    # read rigid body file and generate structure instances for each segment
    listRB=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance)
    
    # score each rigid body segment
    listsc_sccc = []
    print('calculating scores')
    for RB in listRB:
      # sccc score
      score_SCCC=scorer.SCCC(emmap,r,sim_sigma_coeff,structure_instance,RB,c_mode=False)
      SCCC_list_structure_instance.append(score_SCCC)
      print ('>>', score_SCCC)
      listsc_sccc.append(score_SCCC)

    self.assertTrue(len(listRB) == 4)
    self.assertTrue(round(listsc_sccc[0],4) - 0.672 < 0.01) 
    self.assertTrue(round(listsc_sccc[1],4) - 0.642 < 0.01) 
    self.assertTrue(round(listsc_sccc[2],4) - 0.704 < 0.01) 
    self.assertTrue(round(listsc_sccc[3],4) - 0.495 < 0.01) 

  def test_chimera_sccc(self):
    pass

#class MyTestCase2(unittest.TestCase):

if __name__ == '__main__':
  unittest.main()
