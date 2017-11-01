# vim: set expandtab shiftwidth=2 softtabstop=2:

# tool_select.py
#
# Various selection routines for deciding what model / map is selected

from chimerax.core.models import Models
from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure
from chimerax.core.commands import select

def select_model_map(toptool):
  ''' Our way of selecting the model and map for scoring.
  We take the first model we find and the first map we 
  find that are selected. Probably needs improvement.'''

  atomic_model = None
  map_model = None

  for mm in toptool.session.selection.models():
    if isinstance(mm, AtomicStructure):
      atomic_model = mm
      break
  
  for mm in toptool.session.selection.models():
    if isinstance(mm, Volume):
      map_model = mm
      break

  if atomic_model == None or map_model == None:
    toptool._show_error("Please select one model and one map.")
    return (False, None, None)
  
  return (True, atomic_model, map_model)

def select_models_map(toptool):
  ''' Our way of selecting the models and map for scoring with smoc.
  We take the first model we find and the first map we 
  find that are selected. Probably needs improvement.'''

  atomic_models = []
  map_model = None
  
  for mm in toptool.session.selection.models():
    if isinstance(mm, AtomicStructure):
        atomic_models.append(mm)
  
  for mm in toptool.session.selection.models():
    if isinstance(mm, Volume):
        map_model = mm
        break

  if len(atomic_models) == 0 or map_model == None:
    toptool._show_error("Please select one or more models and one map.")
    return (False, None, None)
  
  return (True, atomic_models, map_model)

def select_two(toptool):
  ''' Take the first two things selected. We will see what they are used
  for later. Used in the NMI scoring.'''

  if len(toptool.session.selection.models()) == 2:
    return (True, toptool.session.selection.models()[0], toptool.session.selection.models()[1])

  toptool._show_error("Please select two and only two maps/models.")
  return(False,None,None)


