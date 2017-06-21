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

    # UI consists of a chain selector and search button on top
    # and HTML widget below for displaying results.
    # Layout all the widgets
    from PyQt5.QtWidgets import QLineEdit, QGridLayout, QLabel, QComboBox, QPushButton
    from chimerax.core.ui.widgets import HtmlView
    
    layout = QGridLayout()
    self.top_layout = layout
    layout.setContentsMargins(0, 0, 0, 0)
    # Scoring 
    label = QLabel("SCCC Score:")
    layout.addWidget(label, 0, 0)
    button = QPushButton("calculate")
    button.clicked.connect(self._sccc_score)
    layout.addWidget(button, 0, 2)
    
    # Rigid file
    button_file = QPushButton("select rigid file")
    label_file = QLabel("Rigid filename:")
    self._widget_rigid_file = QLineEdit()
    # TODO - remove this eventually
    self._widget_rigid_file.setText('/home/oni/Projects/ChimeraXTempy/test/rigid_RF.txt')

    button_file.clicked.connect(self._select_rigid_file)
    layout.addWidget(button_file,1,2)
    layout.addWidget(self._widget_rigid_file,1,1)
    layout.addWidget(label_file,1,0)

    label_smoc = QLabel("SMOC Score:")
    layout.addWidget(label_smoc, 2, 0)
    button_smoc = QPushButton("calculate")
    button_smoc.clicked.connect(self._smoc_score)
    layout.addWidget(button_smoc, 2, 2)
    
    #self.results_view = HtmlView(parent, size_hint=(575, 300),
    #                             interceptor=self._html_view,
    #                             schemes=[])
    #layout.addWidget(self.results_view, 1, 0, 1, 3)
    layout.setColumnStretch(0, 0)
    layout.setColumnStretch(1, 10)
    layout.setColumnStretch(2, 0)
    layout.setRowStretch(0, 0)
    layout.setRowStretch(1, 10)
    parent.setLayout(layout)

    # Register for model addition/removal so we can update chain list
    #from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
    #t = session.triggers
    #self._add_handler = t.add_handler(ADD_MODELS, self._update_chains)
    #self._remove_handler = t.add_handler(REMOVE_MODELS, self._update_chains)

  def _sccc_score(self):
    from .sccc import score

    # TODO - Proper checks for existance
    rb_file = self._widget_rigid_file.text()
   
    # Looking for currently selected maps and similar
    print(self.session.models.list())
  
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

    from PyQt5.QtWidgets import QVBoxLayout
    # a figure instance to plot on
    figure = plt.figure()

    # this is the Canvas Widget that displays the `figure`
    # it takes the `figure` instance as a parameter to __init__
    canvas = FigureCanvas(figure)

    parent = self.tool_window.ui_area
    # this is the Navigation widget
    # it takes the Canvas widget and a parent
    toolbar = NavigationToolbar(canvas, parent)
    
    # set the layout
    sublayout = QVBoxLayout()
    sublayout.addWidget(toolbar)
    sublayout.addWidget(canvas)
    
    self.top_layout.addLayout(sublayout,3,0)
    self.top_layout.setColumnStretch(3, 10)

    import random
    ''' plot some random stuff '''
    # random data
    data = [random.random() for i in range(10)]

    # create an axis
    ax = figure.add_subplot(111)

    # discards the old graph
    ax.hold(False)

    # plot data
    ax.plot(data, '*-')

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
