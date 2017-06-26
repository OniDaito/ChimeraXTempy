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

  def _sccc_score(self):
    ''' Run the sccc score as a graphical function, 
    setting the colours of the chosen model.'''

    from .sccc import score

    # TODO - Proper checks for existance
    rb_file = self._widget_rigid_file.text()
   
    # TODO - For now find the first AtomicStructure and first Volume (eventually, do selected and throw errors)
    
    for mm in self.session.models.list():
      if isinstance(mm, AtomicStructure):
        atomic_model = mm
        break
    
    for mm in self.session.models.list():
      if isinstance(mm, Volume):
        map_model = mm
        break
 
    score(self.session, atomic_model, map_model, rb_file)

  def _smoc_score(self):
    ''' Compute the smoc score but also plot
    the scores below the tool.'''

    from .smoc import score
    from PyQt5.QtWidgets import QVBoxLayout

    # TODO - Proper checks for existence
    rb_file = self._widget_rigid_file.text()
  
    # TODO - For now find the first AtomicStructure and first Volume (eventually, do selected and throw errors)
    
    for mm in self.session.models.list():
      if isinstance(mm, AtomicStructure):
        atomic_model = mm
        break
    
    for mm in self.session.models.list():
      if isinstance(mm, Volume):
        map_model = mm
        break

    dict_chains_scores, dict_reslist = score(self.session, atomic_model, map_model, rb_file)

    # a figure instance to plot on
    figure = plt.figure()
    canvas = FigureCanvas(figure)
    parent = self.tool_window.ui_area
    toolbar = NavigationToolbar(canvas, parent)
     
    # TODO - This adds a lot more layers if we keep scoring. A good idea to show improvement
    # perhaps but we may need a way to remove graphs
    sublayout = QVBoxLayout()
    sublayout.addWidget(toolbar)
    sublayout.addWidget(canvas)
    
    self.top_layout.addLayout(sublayout)
  
    # TODO - Really there is only one ch (as one model) for now but eventually there will be more  
    # TODO - multiple models for this score eventually
  
    for ch in dict_chains_scores:
      axes =figure.gca()
      axes.set_ylim([0.4,1.0])
      figure.xlabel = 'Residue_num'
      figure.ylabel = 'SMOC'
      reslist = []
      scorelist = []
      
      for res in dict_reslist[ch]:
        reslist.append(res)
        scorelist.append(dict_chains_scores[ch][res])
     
      ax = figure.add_subplot(111)
      ax.hold(False)
      ax.plot(reslist,scorelist,linewidth=3.0,label="smoc score")
    
    # refresh canvas
    canvas.draw()

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
