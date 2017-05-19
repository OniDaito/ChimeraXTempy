'''
USE THIS SCRIPT TO GENERATE LOCAL CROSS CORRELATION SCORES
E.g. : Initial and output models from FLEX_EM refinement can be scored against the target map.
OUTPUTS local scores and PDB with B-factor records replaced with SMOC scores (CAN BE COLORED IN CHIMERA)
NOTE: that for models with both protein and nucleic acids, the range of local cross correlation values differs for protein and nucleic acid contents and their scores need to be analysed separately.
*NOTE: Residues whose smoc score calculation fails gets a score of 0.0
NOTE: For multiple model inputs, make sure that they have same chain IDs for comparison
'''

import sys
import re
import os
from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.ScoringFunctions import ScoringFunctions
import numpy as np
import shutil,glob
from TEMPy.class_arg import TempyParser

#from datetime import datetime
EXAMPLEDIR = 'Test_Files'
'''
print 'use --help for help'
print '-m/-m1 [map] for input map'
print '-p/-p1 [pdb] for input pdb'
print '-r [resolution]'
print '-t [density threhold], map density threshold (default: 1.5*sigma)'
print '-w [residue window length], length of overlapping residue window used for scoring'
print '--sf [sigma factor], determines the width of the Gaussian used to describe each atom (default 0.187; other values 0.225,0.356,0.425)'
print '\n\n###################################'
'''
#labels for the models [if there are multiple models], fill this if you need customized labels
list_labels = []

tp = TempyParser()
tp.generate_args()
#COMMAND LINE OPTIONS
m1 = tp.args.inp_map1
m2 = tp.args.inp_map2
m = tp.args.inp_map
r1 = tp.args.res1
r2 = tp.args.res2
r = tp.args.res
c1 = tp.args.thr1
c2 = tp.args.thr2
c = tp.args.thr
p = tp.args.pdb
p1 = tp.args.pdb1
p2 = tp.args.pdb2
apix = tp.args.apix
pdir = tp.args.pdir
rigid_out = tp.args.rigidout
rigid_out_prefix = 'rigid'
# the sigma factor determines the width of the Gaussian distribution used to describe each atom
sim_sigma_coeff = 0.187
if not tp.args.sigfac is None:
    sim_sigma_coeff = tp.args.sigfac
#score window
win=9
if not tp.args.window is None:
    win = tp.args.window
#EXAMPLE RUN
flag_example = False
if len(sys.argv) == 1:
  path_example=os.path.join(os.getcwd(),EXAMPLEDIR)
  if os.path.exists(path_example)==True:
    print "%s exists" %path_example
  #else: sys.exit('No input')
  #os.chdir(path_out)
  flag_example = True
  
#GET INPUT DATA
if flag_example:
    map_file = os.path.join(path_example,'1akeA_10A.mrc')
    #p = os.path.join(path_example,'pdb5fuu.ent')
    res_map = 10.0
    DATADIR = path_example
    #MAPDIR = os.path.dir(map_file)
    list_to_check = ['1ake_mdl1.pdb']
    #PDBDIR = DATADIR#os.path.split(os.path.abspath(list_to_check[0]))[0]
elif None in [m1,m2]:
    # for one map and model
    map_file = tp.args.inp_map
    #MAPDIR = os.path.dirname(map_file)
    assert os.path.isfile(map_file)
    print 'reading map'
    Name1 = os.path.basename(m).split('.')[0]
    emmap1=MapParser.readMRC(m)
    if r1 is None and r is None: sys.exit('Input a map and map resolution (required)')
    elif r1 is None: res_map = r
    else: res_map = r1
    if not pdir is None:
        if not os.path.isdir(pdir): sys.exit('Input a single model or directory with multiple conformations')
        DATADIR = pdir
        list_to_check = glob.glob1(DATADIR,'*.pdb')
        #PDBDIR = DATADIR
    elif all(x is None for x in [p,p1,p2]): sys.exit('Input a single model or directory with multiple conformations')
    elif None in [p1,p2]: 
        p = tp.args.pdb
        list_to_check = [os.path.basename(os.path.abspath(p))]
        DATADIR = os.path.dirname(os.path.abspath(p))
        #PDBDIR = os.path.dir(os.path.abspath(p))
    else: sys.exit('Input a single model or directory with multiple conformations')
    

#change these in case of multiple fitted models with single chain
#labels for the models
if len(list_labels) == 0: list_labels = [x.split('.')[0] for x in list_to_check]#['initial','final']
list_styles = [':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-',':','-.','--','-','-']#'--'
#index of model to check z-scores
z_score_check = 2
#----------------------------


def model_tree(list_coord1,distpot=3.5,list_coord2=None):
    try: 
        from scipy.spatial import cKDTree
        coordtree = cKDTree(list_coord2)
    except ImportError:
        from scipy.spatial import KDTree
        coordtree = KDTree(list_coord12)
        #if list_coord2 != None: coordtree1 = KDTree(list_coord2)
    if list_coord2 != None: 
        neigh_points = coordtree.query_ball_point(list_coord1,distpot)
            # use count_neighbors if the corresponding indices are not required
        #else: 
        #    neigh_points = coordtree.query_ball_point(coordtree,distpot)
    #print len(list_coord1), len(neigh_points)
    return neigh_points


start_pdb = list_to_check[0]
iter_num = len(list_to_check)
intermed_file = ""
slow = 0.50
shigh = 0.25 # fraction of structure fitted reasonably well initially
#rigid body file
rigidbody_file = None#

blurrer = StructureBlurrer()
sc = ScoringFunctions()
#read map file
emmap=MapParser.readMRC(map_file)

#-----------------------------
#set plotting parameters
flagplot = 1
try: import matplotlib
except ImportError: flatplot = 0
if flagplot == 1:
  print 'Setting maptpltlib parameters'
  try:
    ##matplotlib.use('Agg')
    try: import matplotlib.pyplot as plt
    except ImportError:flagread = 0 
    from matplotlib import pylab
    try: plt.style.use('ggplot')
    except AttributeError: pass
    legendlist=['iter_0']
    ###colorlist = ['darkcyan','brown','gray','darkgreen','darkkhaki','coral','indigo','maroon']#,'darkkhaki'
    colormap = plt.cm.Spectral#Set1,Spectral#YlOrRd,Spectral,BuGn,Set1,Accent,spring
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 1, iter_num+1)])
    axes = plt.gca()
    plt.gca().set_ylim([0.0,1.0])
    params = {'legend.fontsize': 18}#,'legend.markerscale': 5}
    plt.rcParams['figure.figsize'] = [12, 5]
    plt.rcParams.update(params)
    plt.rcParams.update({'font.size': 18})
    plt.xlabel('Residue Num', fontsize=15)
    plt.ylabel('SMOC',fontsize=15)
    #plt.rcParams['figure.figsize'] = 10, 5
    #plt.figure(figsize=(12,5))
  except: flagplot = 0
#------------------------------
#get scores for start pdb
rfilepath = rigidbody_file
dict_str_scores = {}
if rigidbody_file is not None: rfilepath = os.path.join(DATADIR,rigidbody_file)
#--------------------------------
list_zscores = []
curdir = os.getcwd()
rerun_ct=0
flag_rerun = 0
it = 0
dict_reslist = {}
#flex_em run iter
while iter_num > 0:
  #if it > 0:
  #  prevratio = avghigh/avglow
  #if it > 1: score_inc_prev = score_inc[:]
  ##os.chdir(DATADIR)
  dict_chains_scores = {}
  out_iter_pdb = list_to_check[it]
  lab  = list_labels[it]
  if os.path.isfile(os.path.join(DATADIR,out_iter_pdb)):
    #read pdb
    structure_instance=PDBParser.read_PDB_file('pdbfile',os.path.join(DATADIR,out_iter_pdb),hetatm=False,water=False)
    
    #get scores
    dict_ch_scores,dict_chain_res = sc.SMOC(emmap,res_map,structure_instance,win,rfilepath,sim_sigma_coeff)
  else:
    print 'PDB file not found:', out_iter_pdb
    
  if rigid_out:
    dict_chain_indices, dict_chain_CA = blurrer.get_coordinates(structure_instance)
    rigidf = open(rigid_out_prefix+'_'+lab,'w')
    
  for ch in dict_ch_scores:
    flagch = 1
    dict_res_scores = dict_ch_scores[ch]
    #get res number list (for ref)
    if it == 0:
      dict_reslist[ch] = dict_chain_res[ch][:]
    try: 
      if len(dict_reslist[ch]) == 0: 
        print 'Chain missing:', out_iter_pdb, ch
        flagch = 0
        continue
    except KeyError: 
      print 'Chain not common:', ch, out_iter_pdb
      flagch = 0
      continue
    try: reslist = dict_reslist[ch]
    except KeyError:
      print 'Chain not common:', ch, out_iter_pdb
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
    
    #print list_sccc
    #print scorelist
    #print reslist
    if rigid_out:
      if len(reslist) != len(scorelist):
        print 'Cannot write rigid file for :', lab, ch
      else:
        #if ch in ['',' ',0] and dict_chain_CA.has_key(0): ch_coord = 0
        ch_coord = ch
        #TODO: sort model number issue
        seqres = dict_chain_CA[0][ch_coord].keys()[:]
        seqres.sort()
        list_ca_coord = [dict_chain_CA[0][ch_coord][x] for x in seqres]
        low_list = list_sccc[:int(len(list_sccc)*slow)]
        list_ca_coord_low = [] 
        list_res_low = []     
        for scr in low_list:
          indi = scorelist.index(scr)
          try: resi = int(reslist[indi])
          except TypeError: continue
          #print scr, indi, reslist[indi]
          if dict_chain_CA[0][ch_coord].has_key(resi):
            try: ca_coord = dict_chain_CA[0][ch_coord][resi]
            except KeyError: print 'Cannot write rigid file for :', lab, ch_coord
            list_ca_coord_low.append(ca_coord)
            list_res_low.append(resi)
          else: print resi, dict_chain_CA[0][ch_coord].keys()
        #print list_res_low
        neigh_points = model_tree(list_ca_coord_low,5.0,list_ca_coord)
        list_neigh = []
        for indi in neigh_points:
          for ii in indi:
            list_neigh.append(seqres[ii])
        setlist = set(list_neigh)
        list_flex = list(setlist)
        list_flex.sort()
        #print seqres
        #print list_flex
        for l in range(len(list_flex)):
          try:
            if seqres.index(list_flex[l+1]) - seqres.index(list_flex[l]) != 1:
              rigidf.write(`seqres[seqres.index(list_flex[l])+1]`+':'+ch_coord + ' ' + `seqres[seqres.index(list_flex[l+1])-1]`+':'+ch_coord+"\n")
          except IndexError: pass
        if seqres.index(seqres[-1]) - seqres.index(list_flex[-1]) > 1:
          rigidf.write(`seqres[seqres.index(list_flex[-1])+1]`+':'+ch_coord + ' ' + `seqres[-1]`+':'+ch_coord+"\n")
          
        
    
    '''
    #modify rigid bodies if
    if abs(avghigh/avglow) > 1.0:
      goodset = []
      badset = []
      intrmset = []
      for l in range(len(scorelist)):
        # get indices of good and bad sets based on score
        if scorelist[l] > avghigh1: goodset.append(l)
        elif scorelist[l] < avglow: badset.append(l)
        else: intrmset.append(l)
    else:
      goodset = range(len(scorelist))
    '''
    #print it, 'Num of good scoring residues', len(goodset)
    print list_to_check[it],ch, 'avg-top25%, avg-low25%, avg-high/avg-low', avghigh, avglow, avghigh/avglow
    print list_to_check[it],ch, 'avg', sum(scorelist)/len(scorelist)

    if len(list_to_check) > 0 and len(dict_reslist) == 1 and flagplot == 1:
      plt.plot(reslist,scorelist,linewidth=3.0,label=list_labels[it],linestyle=list_styles[it])

  #include smoc scores as b-factor records
  for x in structure_instance.atomList:
    cur_chain = x.chain
    cur_res = x.get_res_no()
    if not dict_reslist.has_key(cur_chain): continue
    if dict_chains_scores.has_key(cur_chain):
      try: x.temp_fac = dict_chains_scores[cur_chain][cur_res][it]
      except KeyError, IndexError: 
        print 'Residue missing: ',cur_res, ch, out_iter_pdb 
        x.temp_fac = 0.0
    else:
      x.temp_fac = 0.0
    
  #write out pdb file with scores as bfactor records
  pdbid = out_iter_pdb.split('.')[0]
  structure_instance.write_to_PDB(os.path.join(DATADIR,pdbid+"_sc.pdb"))
  it = it+1
  iter_num = iter_num-1
  if rigid_out: rigidf.close()
    #legendlist.append('iter_'+`it`)
#-------------------------------------
##plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
#set legend parameters and plot
if len(list_to_check) > 0 and len(dict_reslist) == 1 and flagplot == 1:
  leg = plt.legend(loc='lower right')
  for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
  plt.show()
elif flagplot == 1: plt.close()

#-------------------------------------
'''
plt.xlabel('Residue Num', fontsize=15)
plt.ylabel('Z-score',fontsize=15)
itns = [z_score_check]
for itn in itns:
  plt.plot(reslist,list_zscores[itn-1],linewidth=3.0,label=list_labels[itn],linestyle=list_styles[itn],color=colorlist[itn])
leg = plt.legend(loc='lower right')
for legobj in leg.legendHandles:
  legobj.set_linewidth(2.0)
plt.show()
'''
try: plt.style.use('ggplot')
except AttributeError: pass
#print dict_str_scores.keys()
for ch in dict_chains_scores:
  it = 0
  axes = plt.gca()
  axes.set_ylim([0.4,1.0])
  plt.xlabel = 'Residue_num'
  plt.ylabel = 'SMOC'
  for model in dict_str_scores:
    if not dict_str_scores[model].has_key(ch): continue
    ch1 = ch
    if ch == ' ': ch1 = ''
    labelname = model+'_'+ch1
    scoreout = os.path.join(DATADIR,"smoc_"+ch1)
    sco = open(scoreout,'w')
    sco.write("Res_num\t"+'\t'.join(list_to_check)+"\n")
    reslist = []
    scorelist = []
    for res in dict_reslist[ch]:
      reslist.append(res)
      scorelist.append(dict_str_scores[model][ch][res])
      sco.write(str(res)+"\t"+'\t'.join(dict_str_scores[model][ch][res])+"\n")
    sco.close() 
    plt.plot(reslist,scorelist,linewidth=3.0,label=labelname,linestyle=list_styles[it],\
             )
    it += 1
  leg = plt.legend(loc='lower right')
  for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
  pylab.savefig(labelname+"_scoreplot.png")
  plt.close()
  
  
