# vim: set expandtab shiftwidth=2 softtabstop=2:

# test_cmd 
#   based on the python3 unittest boilerplate code
#   https://docs.python.org/3.5/library/test.html
#

import unittest

# This might be a bit evil but we need it for now
import sys, os, re
sys.path.insert(0, "../TEMPy")
 
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ShowPlot import Plot
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser

import numpy as np
import shutil,glob

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

class TestSMOC(unittest.TestCase):
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_tempy_smoc(self):
    ''' Test the tempy smoc score based on the files
    provided. Use this as a baseline for the second
    chimeraX test.  It is taken straight from the 
    score_smoc.py example tutorial.'''
  
    list_labels = []

    tp = TempyParser()
    tp.generate_args()
    
    # the sigma factor determines the width of the Gaussian distribution used to describe each atom
    sim_sigma_coeff = 0.187
    #score window
    win=9
     
    path_test = os.getcwd()
    map_file = os.path.join(path_test,'1ake_molmap45_pad10_cubic.mrc')
    res_map = 10.0
    DATADIR = path_test
    list_to_check = ['mdl1.pdb','mdl_ss2.pdb','mdl_ss2_all1.pdb','1ake.pdb']
   
    if len(list_labels) == 0: list_labels = [x.split('.')[0] for x in list_to_check]#['initial','final']
    list_styles = [':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-']#'--'
    
    z_score_check = 2

    def model_tree(list_coord1,distpot=3.5,list_coord2=None):
      try: 
        from scipy.spatial import cKDTree
        coordtree = cKDTree(list_coord2)
      except ImportError:
        from scipy.spatial import KDTree
        coordtree = KDTree(list_coord12)
      if list_coord2 != None: 
        neigh_points = coordtree.query_ball_point(list_coord1,distpot)
      
      return neigh_points

    start_pdb = list_to_check[0]
    iter_num = len(list_to_check)
    intermed_file = ""
    slow = 0.50
    shigh = 0.25 # fraction of structure fitted reasonably well initially
    rigidbody_file = None

    blurrer = StructureBlurrer()
    sc = ScoringFunctions()
    emmap=MapParser.readMRC(map_file)

    rfilepath = rigidbody_file
    dict_str_scores = {}
    if rigidbody_file is not None: rfilepath = os.path.join(DATADIR,rigidbody_file)
    list_zscores = []
    curdir = os.getcwd()
    rerun_ct=0
    flag_rerun = 0
    it = 0
    dict_reslist = {}
   
    # TODO - this whole bit needs a cleanup I think

    while iter_num > 0:
         
      dict_chains_scores = {}
      out_iter_pdb = list_to_check[it]
      lab  = list_labels[it]
      if os.path.isfile(os.path.join(DATADIR,out_iter_pdb)):
        #read pdb
        structure_instance=PDBParser.read_PDB_file('pdbfile',os.path.join(DATADIR,out_iter_pdb),hetatm=False,water=False)
        
        #get scores
        dict_ch_scores,dict_chain_res = sc.SMOC(emmap,res_map,structure_instance,win,rfilepath,sim_sigma_coeff)
      else:
        print('PDB file not found:', out_iter_pdb)
                
      for ch in dict_ch_scores:
        flagch = 1
        dict_res_scores = dict_ch_scores[ch]
        #get res number list (for ref)
        if it == 0:
          dict_reslist[ch] = dict_chain_res[ch][:]
        try: 
          if len(dict_reslist[ch]) == 0: 
            print('Chain missing:', out_iter_pdb, ch)
            flagch = 0
            continue
        except KeyError: 
          print('Chain not common:', ch, out_iter_pdb)
          flagch = 0
          continue
        try: reslist = dict_reslist[ch]
        except KeyError:
          print('Chain not common:', ch, out_iter_pdb)
          flagch = 0
          continue
        if not ch in dict_chains_scores: dict_chains_scores[ch] = {}
        scorelist = []
        for res in reslist:
          try: scorelist.append(dict_res_scores[res])
          except KeyError: 
            if reslist.index(res) <= 0: scorelist.append(dict_res_scores[reslist[reslist.index(res)+1]])
            else: 
              try:  scorelist.append(dict_res_scores[reslist[reslist.index(res)-1]])
              except IndexError: scorelist.append(0.0)
          #save scores for each chain
          curscore = "{0:.2f}".format(round(scorelist[-1],2))
          try: 
            dict_chains_scores[ch][res][it] = str(curscore)
          except KeyError: 
            dict_chains_scores[ch][res] = [str(0.0)]*len(list_to_check)
            dict_chains_scores[ch][res][it] = str(curscore)
        
        dict_str_scores[lab] = dict_chains_scores
            
        #calc ratio between current and prev scores
        if it > 0:
          score_cur = scorelist[:]
          score_inc = [(1+x)/(1+y) for x, y in zip(score_cur, score_prev)][:]
          score_diff = [(x-y) for x, y in zip(score_cur, score_prev)][:]
        #calculate z-scores
        npscorelist = np.array(scorelist)
        try: list_zscores.append((npscorelist-np.mean(npscorelist))/np.std(npscorelist))
        except: list_zscores.append((npscorelist-np.mean(npscorelist)))
        #calculate low and high score bounds
        list_sccc = scorelist[:]
        score_prev = scorelist[:]
        list_sccc.sort()

        #save avg of highest and lowest 20%  
        avglow = list_sccc[int(len(list_sccc)*slow)]
        if avglow == 0.0: avglow = 0.00001
        avghigh = list_sccc[int(len(list_sccc)*(1-shigh))]
        if it == 0: avghigh1 = list_sccc[int(len(list_sccc)*(1-shigh))]
        curratio = avghigh/avglow
                    
        #print it, 'Num of good scoring residues', len(goodset)
        print (list_to_check[it],ch, 'avg-top25%, avg-low25%, avg-high/avg-low', avghigh, avglow, avghigh/avglow)
        print (list_to_check[it],ch, 'avg', sum(scorelist)/len(scorelist))

      #include smoc scores as b-factor records
      for x in structure_instance.atomList:
        cur_chain = x.chain
        cur_res = x.get_res_no()
        if not cur_chain in dict_reslist.keys(): continue
        if cur_chain in dict_chains_scores.keys():
          try: x.temp_fac = dict_chains_scores[cur_chain][cur_res][it]
          except: 
            print ('Residue missing: ',cur_res, ch, out_iter_pdb) 
            x.temp_fac = 0.0
        else:
          x.temp_fac = 0.0
      
      it = it+1
      iter_num = iter_num-1

if __name__ == '__main__':
  unittest.main()
