# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg, FileNameArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf, ModelArg


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


def smoc(session, scoringModel=None, scoringMap=None, rigidFile="rigid.txt", simSigma=0.187, rez=10.0, window=9):
  ''' Calculate the smoc score using the parameters above from the command line
  interface. e.g smoc (#1) (#2) rigid.txt . Optional parameters are the sigma, rez and 
  window.'''
  
  from .smoc import score
  score(session, scoringModel, scoringMap, rigidFile,simSigma,rez,window)

  
smoc_desc = CmdDesc(required=[("scoringModel", ModelArg), 
                          ("scoringMap", ModelArg),
                          ("rigidFile", FileNameArg),],
                          optional=[("simSigma",FloatArg),
                            ("rez",FloatArg),
                            ("window",IntArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SMOC function")
    
def register_smoc():
    from chimerax.core.commands import register
    register('smoc', smoc_desc, smoc)
