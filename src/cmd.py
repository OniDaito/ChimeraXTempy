# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg, OpenFileNameArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf, ModelArg, ModelsArg
from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

#
# SCCC Section
#
def sccc(session, scoringModel=None, scoringMap=None, rigidFile="none", rez=10.0, simSigma=0.187 ):
  ''' Calculate the sccc score using the parameters above from the command line
  interface. e.g sccc (#1) (#2) rigid.txt 10.0. Optional parameters are the sigma.'''
  
  from .sccc import score
  score(session, scoringModel, scoringMap, rigidFile, rez, simSigma)
  
sccc_desc = CmdDesc(required=[("scoringModel", ModelArg), 
                          ("scoringMap", ModelArg),
                          ("rigidFile", OpenFileNameArg),("rez",FloatArg)],
                           optional=[("simSigma",FloatArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SCCC function")
    
def register_sccc():
    from chimerax.core.commands import register
    register('sccc', sccc_desc, sccc)

#
# CCC Section
#
def ccc(session, scoringModel=None, scoringMap=None, rez=10.0):
  ''' Calculate the sccc score using the parameters above from the command line
  interface. e.g ccc (#1) (#2) 10.0'''
  
  from .ccc import score
  score(session, scoringModel, scoringMap, rez)

ccc_desc = CmdDesc(required=[("scoringModel", ModelArg), 
                          ("scoringMap", ModelArg),
                          ("rez",FloatArg)],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy CCC function")
    
def register_ccc():
    from chimerax.core.commands import register
    register('ccc', ccc_desc, ccc)


#
# SMOC 
#
def smoc(session, scoringModels=None, scoringMap=None, rez=10.0, simSigma=0.187, window=9):
  ''' Calculate the smoc score using the parameters above from the command line
  interface. e.g smoc (#1 #2) (#3). Optional parameters are the sigma and 
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
  scores = score(session, scoringModels, scoringMap, "", rez, simSigma,window)

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
                          ("scoringMap", ModelArg),("rez",FloatArg),],
                          optional=[("simSigma",FloatArg),
                            ("window",IntArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SMOC function")
    
def register_smoc():
    from chimerax.core.commands import register
    register('smoc', smoc_desc, smoc)

#
# NMI for single items
#
def nmi(session, scoringMapModel1=None, scoringMapModel2=None, rez1=None, rez2=None, contour1=None, contour2=None):
  ''' Calculate the nmi score using the parameters above from the command line
  interface. e.g nmi (#1) (#2) 6.6 4.0. Optional parameters are the resolutions and contours.'''
  
  from .nmi import score
  from PyQt5 import QtWidgets

  nmi_score = 0.0
  if isinstance(scoringMapModel1, Volume) and isinstance(scoringMapModel2, AtomicStructure):
    nmi_score = score(session, scoringMapModel2, scoringMapModel1, None, None, rez1, rez2, contour1, contour2 )
  elif isinstance(scoringMapModel1, Volume) and isinstance(scoringMapModel2, Volume):
    nmi_score = score(session, None, scoringMapModel1, None, scoringMapModel2, rez1, rez2, contour1, contour2 )
  else :
    print("TEMPY Error: Please provide a map and model or two maps.")
    return
    
  print("NMI Score: ", nmi_score)

#nmi_desc = CmdDesc(required=[("scoringMapModel1", ModelArg), 
#                          ("scoringMapModel2", ModelArg),],
#                          optional=[("rez1",FloatArg),
#                            ("rez2",FloatArg),
#                            ("contour1",FloatArg),
#                            ("contour2",FloatArg),],
#                      keyword=[("log", BoolArg)],
#                      synopsis="The Tempy NMI function")
    
#def register_nmi():
#    from chimerax.core.commands import register
#    register('nmi', nmi_desc, nmi)


#
# difmap for two maps
#
def difmap(session, map0 = None, map1 = None, rez1 = None, rez2 = None, contour1 = None, contour2 = None ):
  ''' Calculate the nmi score using the parameters above from the command line
  interface. e.g difmap (#1) (#2) 10 10'''
  
  from .difmap import genmap
  from PyQt5 import QtWidgets

  genmap(session, map0, map1, rez1, rez2, contour1, contour2)

difmap_desc = CmdDesc(required=[("map0", ModelArg), 
                          ("map1", ModelArg),("rez1",FloatArg),
                            ("rez2",FloatArg),],
                          optional=[
                            ("contour1",FloatArg),
                            ("contour2",FloatArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy difmap function")
    
def register_difmap():
    from chimerax.core.commands import register
    register('difmap', difmap_desc, difmap)

#
# NMI for multiple objects
#
def nmi(session, comparators = None, compared = None, rez1 = None, rez2 = None, contour1 = None, contour2 = None ):
  ''' Calculate the nmi score using the parameters above from the command line
  interface. e.g nmi (#1 #2) (#3) 6.6 4.0 10.0 10.0. Optional parameters are the contours.
  We assume ALL comparators have the same rez and contours.'''
  
  from .nmi import score_cmd
  from PyQt5 import QtWidgets

  nmi_scores = []
  contour_compared = contour2
  rez_compared = rez2
  contours_comparators = []
  rez_comparators = []

  for c in comparators:
    rez_comparators.append(rez1)
    contours_comparators.append(contour1)

  nmi_scores = score_cmd(session, comparators, compared, rez_comparators, rez_compared, contours_comparators, contour_compared)

  if nmi_scores != None: 
    for nmi_score in nmi_scores:
      print("NMI Score: ", nmi_score)

nmi_desc = CmdDesc(required=[("comparators", ModelsArg), 
                          ("compared", ModelArg),("rez1",FloatArg),
                            ("rez2",FloatArg),],
                          optional=[
                            ("contour1",FloatArg),
                            ("contour2",FloatArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy NMI function")
    
def register_nmi():
    from chimerax.core.commands import register
    register('nmi', nmi_desc, nmi)


