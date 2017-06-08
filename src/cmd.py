# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf


def sccc_score(session, spec=None, log=None):
  print("Calculate SCCC Score")
  
  
sccc_score_desc = CmdDesc(required=[("spec", AtomSpecArg),],
                        keyword=[("log", BoolArg)],
                        synopsis="The Tempy Program")
    

# I suspect this refers to the cmd line bit
def register_sccc_score():
    from chimerax.core.commands import register, CmdDesc, NoArg
    desc = CmdDesc(
        required  = [],
        optional  = [('sim', NoArg),] 
                )
    register('sccc_score', desc, sccc_score)
