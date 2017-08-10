# vim: set expandtab shiftwidth=2 softtabstop=2:

# test_cmd 
#   based on the python3 unittest boilerplate code
#   https://docs.python.org/3.5/library/test.html
#

import unittest

# This might be a bit evil but we need it for now
import sys, os, re
sys.path.insert(0, "../TEMPy")
 
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ShowPlot import Plot
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser
from TEMPy.StructureBlurrer import StructureBlurrer

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
    p = os.path.join(path_test,'1ake_mdl1.pdb')
    r = 10.0
    rb_file = os.path.join(path_test,'1ake_mdl1_rigid.txt')
    
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
    
    for RB in listRB:
      # sccc score
      score_SCCC=scorer.SCCC(emmap,r,sim_sigma_coeff,structure_instance,RB)
      listsc_sccc.append(score_SCCC)

    self.assertTrue(len(listRB) == 6)
    self.assertTrue(abs(round(listsc_sccc[0],4) - 0.954) < 0.01) 
    self.assertTrue(abs(round(listsc_sccc[1],4) - 0.427) < 0.01) 
    self.assertTrue(abs(round(listsc_sccc[2],4) - 0.624) < 0.01) 
    self.assertTrue(abs(round(listsc_sccc[3],4) - 0.838) < 0.01) 
    self.assertTrue(abs(round(listsc_sccc[4],4) - 0.971) < 0.01) 
    self.assertTrue(abs(round(listsc_sccc[5],4) - 0.928) < 0.01) 

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
    map_file = os.path.join(path_test,'1akeA_10A.mrc')
    res_map = 10.0
    DATADIR = path_test
    list_to_check = ['1ake_mdl1.pdb']
   
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
                    
        self.assertTrue( abs(avghigh - 0.967) < 0.01)
        self.assertTrue( abs(avglow - 0.956) < 0.01)
        self.assertTrue( abs(sum(scorelist)/len(scorelist) - 0.899) < 0.01)


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


#calculate map contour
def map_contour(m,t=-1.):
  mName = os.path.basename(m).split('.')[0]
  #print 'reading map'
  emmap=MapParser.readMRC(m)
  c1 = None
  if t != -1.0:
    zeropeak,ave,sigma1 = emmap._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return mName,emmap,c1

#calculate model contour
def model_contour(p,res=4.0,emmap=False,t=-1.):
  pName,modelmap = blur_model(p,res,emmap)
  c1 = None
  if t != -1.0:
    c1 = t*emmap.std()#0.0
  return pName,modelmap,c1

def blur_model(p,res=4.0,emmap=False):
  pName = os.path.basename(p).split('.')[0]
  structure_instance=PDBParser.read_PDB_file(pName,p,hetatm=False,water=False)
  blurrer = StructureBlurrer()
  if res is None:
    sys.exit('Map resolution required..')
  modelmap = blurrer.gaussian_blur_real_space(structure_instance, res,densMap=emmap,normalise=True) 
  return pName,modelmap

def lpfilter(emmap,r):
  cutoff = emmap.apix/float(r)
  mapfilt = emmap._tanh_lowpass(cutoff)
  return mapfilt 

def match_grid(emmap1,emmap2,c1,c2):
  # DETERMINE A COMMON ALIGNMENT BOX : fill minvalue for extra voxel pads
  spacing = emmap2.apix
  if emmap2.apix < emmap1.apix: spacing = emmap1.apix
  grid_shape, new_ori = emmap1._alignment_box(emmap2,spacing)
  # INTERPOLATE TO NEW GRID
  try: emmap_1 = emmap1._interpolate_to_grid1(grid_shape,spacing,new_ori)
  except: emmap_1 = emmap1._interpolate_to_grid(grid_shape,spacing,new_ori)
  try: c1 = emmap_1._find_level(np.sum(emmap1.fullMap>c1)*(emmap1.apix**3))
  except: pass
  del emmap1.fullMap
  del emmap1
  try: emmap_2 = emmap2._interpolate_to_grid1(grid_shape,spacing,new_ori)
  except: emmap_2 = emmap2._interpolate_to_grid(grid_shape,spacing,new_ori)
  try: c2 = emmap_2._find_level(np.sum(emmap2.fullMap>c2)*(emmap2.apix**3))
  except: pass
  del emmap2.fullMap
  del emmap2
  return emmap_1, emmap_2



class TestNMI(unittest.TestCase):
  # Only use setUp() and tearDown() if necessary
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_tempy_nmi(self):
    ''' Test the tempy nmi score based on the files
    provided. Use this as a baseline for the second
    chimeraX test. '''

    path_test = "./"
    m = os.path.join(path_test,'emd_5168.map')
    p = os.path.join(path_test,'emd_5170.map')
   
    sc = ScoringFunctions()

    rez1 = 6.6
    rez2 = 15.0

    Name1,emmap1,c1 = map_contour(m,t=1.5)
    Name2,emmap2,c2 = map_contour(p,t=1.5)
   
    print(rez1,rez2,c1,c2,emmap1.apix, emmap2.apix)

    if not sc.mapComparison(emmap1,emmap2):
      emmap1._crop_box(c1,0.5)
      emmap2._crop_box(c2,0.5)
      
      if rez1 > 1.25*rez2: 
        emmap_2 = lpfilter(emmap2,rez1)
        emmap1, emmap2 = match_grid(emmap1,emmap_2,c1,c2)
      elif rez2 > 1.25*rez1:
        emmap_1 = lpfilter(emmap1,rez2)
        emmap1, emmap2 = match_grid(emmap_1,emmap2,c1,c2)
      else:
        emmap1, emmap2 = match_grid(emmap1,emmap2,c1,c2)
   
    nmi = 0
    try:
      nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
      if nmi < 0.0: nmi = 0.0
    except:
      self.assertTrue(False)
      print_exc()
      nmi = 0.0
   
    self.assertTrue(abs(round(nmi,5) - 1.0492) < 0.001) 

    # Now test with a model and map
    p = os.path.join(path_test,'1J6Z.pdb')
    m = os.path.join(path_test,'emd_5168_monomer.mrc')
    res = 6.6
    Name1 = os.path.basename(m).split('.')[0]
    Name2 = os.path.basename(p).split('.')[0]
    emmap1=MapParser.readMRC(m)
    structure_instance=PDBParser.read_PDB_file(Name2,p,hetatm=False,water=False)
    blurrer = StructureBlurrer()
    emmap2 = blurrer.gaussian_blur(structure_instance, res,densMap=emmap1)
    c1 = 9.7
    c2 = 1.0

    nmi = 0
    try:
      nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
      if nmi < 0.0: nmi = 0.0
    except:
      self.assertTrue(False)
      print_exc()
      nmi = 0.0
  
    self.assertTrue(abs(round(nmi,5) - 1.0575) < 0.001) 

if __name__ == '__main__':
  unittest.main()
