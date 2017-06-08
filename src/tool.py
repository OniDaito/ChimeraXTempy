# vim: set expandtab shiftwidth=2 softtabstop=2:

# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
# Seems like Tool exists inside the bundle and needs overridding


from chimerax.core.tools import ToolInstance
from chimerax.core.models import Models

from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

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
