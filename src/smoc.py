# vim: set expandtab shiftwidth=2 softtabstop=2:

# sccc
#   "score" - compute the sccc score using tempy
#

from TEMPy.MapParser import MapParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom

from .util import chimera_to_tempy_atom, chimera_to_tempy_map

import numpy as np

import math, os

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


def score(session, atomic_model, map_model, rigid_filename, colour_atoms=True):
  atomlist = []

  # TODO - rigid_filename might be optional?

  sc = ScoringFunctions()
  
  for atom in atomic_model.atoms:
    atomlist.append(chimera_to_tempy_atom(atom, len(atomlist)))

  bio_atom_structure = BioPy_Structure(atomlist)
  bio_map_structure = chimera_to_tempy_map(map_model)

  sim_sigma_coeff = 0.187
  win=9
  res_map = 10.0
  slow = 0.50
  shigh = 0.25 # fraction of structure fitted reasonably well initially
  list_zscores = []
  curdir = os.getcwd()
  rerun_ct=0
  flag_rerun = 0
  it = 0
  dict_reslist = {}
 
  # TODO - So I've removed the iteration bit because we need to think a little more
  # about what PDBs to check I suspect

  dict_chains_scores = {}


  dict_ch_scores,dict_chain_res = sc.SMOC(bio_map_structure, res_map, bio_atom_structure, win, rigid_filename, sim_sigma_coeff)
 
  # Option to save the rigid file?
  #if rigid_out:
  #  dict_chain_indices, dict_chain_CA = blurrer.get_coordinates(structure_instance)
  #  rigidf = open(rigid_out_prefix+'_'+lab,'w')
    
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
        dict_chains_scores[ch][res] = [str(0.0)]
        dict_chains_scores[ch][res][it] = str(curscore)
    
        
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
    print (ch, 'avg-top25%, avg-low25%, avg-high/avg-low', avghigh, avglow, avghigh/avglow)
    print (ch, 'avg', sum(scorelist)/len(scorelist))

    # if len(list_to_check) > 0 and len(dict_reslist) == 1 and flagplot == 1:
    #   plt.plot(reslist,scorelist,linewidth=3.0,label=list_labels[it],linestyle=list_styles[it])

    print("seem to have finished a bit")

    #include smoc scores as b-factor records
    for x in bio_atom_structure.atomList:
      cur_chain = x.chain
      cur_res = x.get_res_no()
      if not cur_chain in dict_reslist.keys(): continue
      if cur_chain  in dict_chains_scores.keys():
        try: x.temp_fac = dict_chains_scores[cur_chain][cur_res][it]
        except: 
          print ('Residue missing: ',cur_res, ch, out_iter_pdb) 
          x.temp_fac = 0.0
      else:
        x.temp_fac = 0.0
      
  #-------------------------------------
  ##plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
  #set legend parameters and plot
  #if len(list_to_check) > 0 and len(dict_reslist) == 1 and flagplot == 1:
  # leg = plt.legend(loc='lower right')
  #  for legobj in leg.legendHandles:
  #    legobj.set_linewidth(2.0)
  #  plt.show()
  #elif flagplot == 1: plt.close()

  #-------------------------------------
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
    '''


