# vim: set expandtab shiftwidth=2 softtabstop=2:

# This file seems to hold the actual commands we can execute inside
# of ChimeraX. Its called from __init__.py I believe

from chimerax.core.commands import CmdDesc, AtomSpecArg
from chimerax.core.commands import StringArg, BoolArg, FloatArg, IntArg, EnumOf


def tempy(session, spec=None, log=None):
  print("Testing tempy")
  
  
tempy_desc = CmdDesc(required=[("spec", AtomSpecArg),],
                        keyword=[("log", BoolArg)],
                        synopsis="The Tempy Program")
    
    
def register_tempy():
    from chimerax.core.commands import register, CmdDesc, NoArg
    desc = CmdDesc(
        required  = [],
        optional  = [('sim', NoArg),] 
                )
    register('tempy', desc, tempy)


