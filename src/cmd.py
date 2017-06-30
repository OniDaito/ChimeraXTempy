# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg, FileNameArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf, ModelArg, ModelsArg
from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure


def sccc(session, scoringModel=None, scoringMap=None, rigidFile="rigid.txt", simSigma=0.187, rez=10.0):
  ''' Calculate the sccc score using the parameters above from the command line
  interface. e.g sccc (#1) (#2) rigid.txt . Optional parameters are the sigman and rez.'''
  
  from .sccc import score
  score(session, scoringModel, scoringMap, rigidFile, simSigma, rez)

  
sccc_desc = CmdDesc(required=[("scoringModel", ModelArg), 
                          ("scoringMap", ModelArg),
                          ("rigidFile", FileNameArg),],
                           optional=[("simSigma",FloatArg),
                            ("rez",FloatArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SCCC function")
    

def register_sccc():
    from chimerax.core.commands import register
    register('sccc', sccc_desc, sccc)


def smoc(session, scoringModels=None, scoringMap=None, simSigma=0.187, rez=10.0, window=9):
  ''' Calculate the smoc score using the parameters above from the command line
  interface. e.g smoc (#1,#2) (#3). Optional parameters are the sigma, rez and 
  window.'''
  
  from .smoc import score
  from PyQt5 import QtWidgets

  # Check the parameters for models and maps
  for mm in scoringModels:
    if not isinstance(mm, AtomicStructure):
      print("TEMPY Error: Please provide model(s) for the first parameter.")
      return
    
    if not isinstance(scoringMap, Volume):
      print("TEMPY Error: Please provide a map for the second parameter.")
      return
   
  # Now do the scoring and print the scores
  scores = score(session, scoringModels, scoringMap, "",simSigma,rez,window)

  idx = 1
  for  (dict_chains_scores, dict_reslist) in scores:
    # So this is a text based approach so we print out each score
    print("SMOC #" + str(idx))

    for ch in dict_chains_scores:
      reslist = []
      scorelist = []
        
      for res in dict_reslist[ch]:
        reslist.append(res)
        scorelist.append(dict_chains_scores[ch][res])
       
    for i in range(0,len(reslist)):
      print(reslist[i], scorelist[i])
    
    print("---")
    print()
    idx+=1


smoc_desc = CmdDesc(required=[("scoringModels", ModelsArg), 
                          ("scoringMap", ModelArg),],
                          optional=[("simSigma",FloatArg),
                            ("rez",FloatArg),
                            ("window",IntArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SMOC function")
    
def register_smoc():
    from chimerax.core.commands import register
    register('smoc', smoc_desc, smoc)
