# vim: set expandtab shiftwidth=2 softtabstop=2:

# tool_scores.py
#
# The functions for computing the scores when the corresponding
# button in the tool is pressed

from . import tool_select
from . import tool_layout

import matplotlib
matplotlib.use("Qt5Agg", force=True)

from chimerax.core.models import Models
from chimerax.core.commands import select
from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import os

def smoc_score(toptool):
  ''' Check all the variables and perform the smoc score, creating
  the corresponding graph in the master tool.'''
  from .smoc import score
  from PyQt5.QtWidgets import QVBoxLayout

  # Check Rigid file
  # Optional with smoc
  rb_file = toptool._widget_rigid_file_smoc.text()
  if not os.path.isfile(rb_file):
    rb_file =""

  # Check model and map
  result, atomic_models, map_model = tool_select.select_models_map(toptool)    
  if not result:
    return

  # Check the options
  try:
    sim_sigma = float(toptool._widget_sigma_smoc.text())
    rez = float(toptool._widget_rez_smoc.text())
    win = int(float(toptool._widget_window_smoc.text()))
  except:
    toptool._show_error("Check the values for rez, sigma and window.")
    return

  # a figure instance to plot on
  if toptool._figure == None:
    toptool._figure = matplotlib.figure.Figure()

  # TODO - This adds a lot more layers if we keep scoring. A good idea to show improvement
  # perhaps but we may need a way to remove graphs. For now, lets just go with replacement. 
  if toptool._canvas == None:
    toptool._canvas = FigureCanvas(toptool._figure)
    
    parent = toptool.tool_window.ui_area
    toolbar = NavigationToolbar(toptool._canvas, parent)
   
    sublayout = QVBoxLayout()
    sublayout.addWidget(toolbar)
    sublayout.addWidget(toptool._canvas)
  
    toptool.top_layout.addLayout(sublayout)
    toptool._subplot = toptool._figure.add_subplot(111)
    toptool._figure.xlabel = 'Residue_num'
    toptool._figure.ylabel = 'SMOC'
    toptool._canvas.mpl_connect('button_press_event', toptool.onclick)

  toptool._subplot.cla()
  toptool._figure.subplots_adjust(bottom=0.25)
 
  # Call score
  idx = 0
  toptool._current_smoc_scores = [] # keep for graph clicking
  toptool._scored_models_smoc = atomic_models
  
  for (dict_chains_scores, dict_reslist) in score(toptool.session, atomic_models,
      map_model, rb_file, rez, sim_sigma, win):

    toptool._current_smoc_scores.append([])

    # TODO - we are mixing chains here. Need to check that sort of thing really
    for ch in dict_chains_scores:
      reslist = []
      scorelist = []
      
      for res in dict_reslist[ch]:
        reslist.append(res)
        tp = dict_chains_scores[ch][res] 
        scorelist.append(tp)
        toptool._current_smoc_scores[idx].append(tp)

      col = atomic_models[idx].single_color
      col = (float(col[0])/256.0,float(col[1])/256.0,float(col[2])/256.0)
      toptool._subplot.plot(reslist,scorelist,linewidth=1.0,label="smoc score", color=col)
    
    idx+=1

  # refresh canvas
  toptool._canvas.draw()

def sccc_score(toptool):
  ''' Run the sccc score as a graphical function, 
  setting the colours of the chosen model.'''

  from .sccc import score

  # Check rigid score file
  rb_file = toptool._widget_rigid_file_sccc.text()
  if not os.path.isfile(rb_file):
    toptool._show_error("File " + rb_file + " does not exist")
    return

  try:
    sim_sigma = float(toptool._widget_sigma_sccc.text())
    rez = float(toptool._widget_rez_sccc.text())
  except:
    toptool._show_error("Check the values for rez and sigma")
    return

  # Find models
  result, atomic_model, map_model = tool_select.select_model_map(toptool)
  if result:
    scores = score(toptool.session, atomic_model, map_model, rb_file, rez, sim_sigma)
    tool_layout.add_sccc_scores(toptool, scores) 
    # Sadly we need to keep current scores so we can reference them elsewhere :/
    toptool._current_sccc_scores = scores
    toptool._scored_models_sccc = atomic_model
 

def nmi_score(toptool):
  ''' Run the nmi score, printing to the log'''
  from .nmi import score

  # Check contours
  contour1 = 1.5
  contour2 = 1.5
  if toptool._widget_c1_nmi.isEnabled() and toptool._widget_c1_nmi.isEnabled():        
    try:
      contour1 = float(toptool._widget_c1_nmi.text())
      contour2 = float(toptool._widget_c2_nmi.text()) 
    except:
      toptool._show_error("Check the values fo contour1 and contour2")
      return

  # Find models
  result, scoringMapModel1, scoringMapModel2 = tool_select.select_two(toptool)
  nmi_score = 0.0
  if result:
    if isinstance(scoringMapModel1, Volume) and isinstance(scoringMapModel2, AtomicStructure):
      try:
        rez1 = float(toptool._widget_rez1_nmi.text()) 
        rez2 = float(toptool._widget_rez2_nmi.text())
        nmi_score = score(toptool.session, scoringMapModel2, scoringMapModel1, None, None, rez1, rez2, contour1, contour2 )
      except: 
        toptool._show_error("Check the values for rez1 and rez2")
        return

    elif isinstance(scoringMapModel2, Volume) and isinstance(scoringMapModel1, AtomicStructure):
      try:
        rez1 = float(toptool._widget_rez1_nmi.text())
        rez2 = float(toptool._widget_rez2_nmi.text())
        nmi_score = score(toptool.session, scoringMapModel1, scoringMapModel2, None, None, rez1, rez2, contour1, contour2 )
      except:
        toptool._show_error("Check the values for rez1 and rez2")
        return

    elif isinstance(scoringMapModel1, Volume) and isinstance(scoringMapModel2, Volume):
       try:
        rez1 = float(toptool._widget_rez1_nmi.text())
        rez2 = float(toptool._widget_rez2_nmi.text())
        nmi_score = score(toptool.session, None, scoringMapModel1, None, scoringMapModel2, rez1, rez2, contour1, contour2 )
       except:
        toptool._show_error("Check the values for rez1 and rez2")
        return
    else :
      toptool._show_error("Please provide a map and model, or two maps.")
      return
      
    print("NMI Score: ", nmi_score)


