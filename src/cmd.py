# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg, FileNameArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf, ModelArg


def sccc(session, scoringModel=None, scoringMap=None, rigidFile="rigid.txt"):
  ''' Calculate the sccc score using the parameters above from the command line
  interface. e.g sccc (#1) (#2) rigid.txt '''
  
  from .sccc import score
  score(session, scoringModel, scoringMap, rigidFile)

  
sccc_desc = CmdDesc(required=[("scoringModel", ModelArg), 
                          ("scoringMap", ModelArg),
                          ("rigidFile", FileNameArg),],
                      keyword=[("log", BoolArg)],
                      synopsis="The Tempy SCCC function")
    

# I suspect this refers to the cmd line bit
def register_sccc():
    from chimerax.core.commands import register
    register('sccc', sccc_desc, sccc)
