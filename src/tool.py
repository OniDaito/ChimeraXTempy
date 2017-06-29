# vim: set expandtab shiftwidth=2 softtabstop=2:

# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
# Seems like Tool exists inside the bundle and needs overridding

from chimerax.core.tools import ToolInstance
from chimerax.core.models import Models

from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure


from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pyplot as plt
import os

class ToolUI(ToolInstance):

  def __init__(self, session, tool_name):

    # Standard template stuff
    ToolInstance.__init__(self, session, tool_name)
    self.display_name = "Tempy"
  
    from chimerax.core.ui.gui import MainToolWindow
    self.tool_window = MainToolWindow(self)
    self.tool_window.manage(placement="side")
    parent = self.tool_window.ui_area

    # All the PyQT5 Layout business
    from PyQt5.QtWidgets import QLineEdit, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QPushButton
    from chimerax.core.ui.widgets import HtmlView
    
    layout = QVBoxLayout()
    self.top_layout = layout
    layout.setContentsMargins(0, 0, 0, 0)

    # Score layouts
    button_layout = QHBoxLayout()
    layout.addLayout(button_layout)

    # Scoring 
    button_sccc = QPushButton("SCCC")
    button_sccc.clicked.connect(self._sccc_score)
    button_layout.addWidget(button_sccc)
   
    button_smoc = QPushButton("SMOC")
    button_smoc.clicked.connect(self._smoc_score)
    button_layout.addWidget(button_smoc)
  
    # Options for the two scores

    options_layout = QHBoxLayout()
    layout.addLayout(options_layout)

    label_rez = QLabel("Res.")
    label_sigma = QLabel("Sigma")
    label_window = QLabel("Window")

    self._widget_sigma = QLineEdit()
    self._widget_rez = QLineEdit()
    self._widget_window = QLineEdit()

    options_layout.addWidget(label_rez)
    options_layout.addWidget(self._widget_rez)
    options_layout.addWidget(label_sigma)
    options_layout.addWidget(self._widget_sigma)
    options_layout.addWidget(label_window)
    options_layout.addWidget(self._widget_window)

    self._widget_rez.setText("10.0")
    self._widget_sigma.setText("0.187")
    self._widget_window.setText("9")

    # Rigid file
    rigid_layout = QHBoxLayout()
    layout.addLayout(rigid_layout)

    button_file = QPushButton("Rigid file")
    self._widget_rigid_file = QLineEdit()
    
    # TODO - remove this eventually
    self._widget_rigid_file.setText('/home/oni/Projects/ChimeraXTempy/test/rigid_RF.txt')

    button_file.clicked.connect(self._select_rigid_file)
    rigid_layout.addWidget(button_file)
    rigid_layout.addWidget(self._widget_rigid_file)

    parent.setLayout(layout)

    # Figures for SMOC
    self._figure = None
    self._canvas = None

  def _select_model_map(self):
    ''' Our way of selecting the model and map for scoring.
    We take the first model we find and the first map we 
    find that are selected. Probably needs improvement.'''

    atomic_model = None
    map_model = None

    for mm in self.session.models.list():
      if isinstance(mm, AtomicStructure):
        # TODO - selected doesnt seem to work too well :/
        #if mm.selected:
        atomic_model = mm
        break
    
    for mm in self.session.models.list():
      if isinstance(mm, Volume):
        if mm.selected:
          map_model = mm
          break
 
    if atomic_model == None or map_model == None:
      print("TEMPY Error: Please select one model and one map.")
      return (False, None, None)
    
    return (True, atomic_model, map_model)

  def _sccc_score(self):
    ''' Run the sccc score as a graphical function, 
    setting the colours of the chosen model.'''

    from .sccc import score

    # Check rigid score file
    rb_file = self._widget_rigid_file.text()
    if not os.path.isfile(rb_file):
      print("TEMPY error: File " + rb_file + " does not exist")
      return

    try:
      sim_sigma = float(self._widget_sigma.text())
      rez = float(self._widget_rez.text())
    except:
      print("TEMPY Error: Check the values for rez and sigma")
      return

    # Find models
    result, atomic_model, map_model = self._select_model_map()
    if result:
      score(self.session, atomic_model, map_model, rb_file, sim_sigma, rez)

  def _smoc_score(self):
    ''' Compute the smoc score but also plot
    the scores below the tool.'''

    from .smoc import score
    from PyQt5.QtWidgets import QVBoxLayout

    # Check Rigid file
    # Optional with smoc
    rb_file = self._widget_rigid_file.text()
    if not os.path.isfile(rb_file):
      rb_file =""

    # Check model and map
    result, atomic_model, map_model = self._select_model_map()    
    if not result:
      return
  
    # Check the options
    try:
      sim_sigma = float(self._widget_sigma.text())
      rez = float(self._widget_rez.text())
      win = int(float(self._widget_window.text()))
    except:
      print("TEMPY Error: Check the values for rez, sigma and window.")
      return
 
    # Call score
    dict_chains_scores, dict_reslist = score(self.session, atomic_model,
          map_model, rb_file, sim_sigma, rez, win)

    # a figure instance to plot on
    if self._figure == None:
      self._figure = plt.figure()
   
    # TODO - This adds a lot more layers if we keep scoring. A good idea to show improvement
    # perhaps but we may need a way to remove graphs. For now, lets just go with replacement. 
    if self._canvas == None:
      self._canvas = FigureCanvas(self._figure)
      parent = self.tool_window.ui_area
      toolbar = NavigationToolbar(self._canvas, parent)
     
      sublayout = QVBoxLayout()
      sublayout.addWidget(toolbar)
      sublayout.addWidget(self._canvas)
    
      self.top_layout.addLayout(sublayout)
      self._subplot = self._figure.add_subplot(111)


    # TODO - Really there is only one ch (as one model) for now but eventually there will be more. This all depends on how we work out the selection stuff. Maybe use the model window? 
    for ch in dict_chains_scores:
      self._figure.xlabel = 'Residue_num'
      self._figure.ylabel = 'SMOC'
      reslist = []
      scorelist = []
      
      for res in dict_reslist[ch]:
        reslist.append(res)
        scorelist.append(dict_chains_scores[ch][res])
     
      self._subplot.cla()
      self._subplot.hold(False)
      self._subplot.plot(reslist,scorelist,linewidth=3.0,label="smoc score")
    
    # refresh canvas
    self._canvas.draw()

  def _select_rigid_file(self):
    from PyQt5.QtWidgets import QFileDialog
    filename = QFileDialog.getOpenFileName(None, 'OpenFile')
    print(filename)
    self._widget_rigid_file.setText(filename[0])

  def _html_view(self):
    print("HTML Test")

  def delete(self):
    #t = self.session.triggers
    #t.remove_handler(self._add_handler)
    #t.remove_handler(self._remove_handler)
    super().delete()

  def take_snapshot(self, session, flags):
    # For now, do not save anything in session.
    # Need to figure out which attributes (like UI widgets)
    # should start with _ so that they are not saved in sessions.
    # And add addition data to superclass data.
    return super().take_snapshot(session, flags)

  @classmethod
  def restore_snapshot(cls, session, data):
    # For now do nothing.  Should unpack data and restart tool.
    return None
