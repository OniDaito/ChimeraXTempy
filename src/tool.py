# vim: set expandtab shiftwidth=2 softtabstop=2:

# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
# Seems like Tool exists inside the bundle and needs overridding

from chimerax.core.tools import ToolInstance
from PyQt5.QtWidgets import QMessageBox

import os, math 

from . import tool_layout
from . import tool_scores

class ToolUI(ToolInstance):

  def __init__(self, session, tool_name):
    ''' Call the special layout function '''
    ToolInstance.__init__(self, session, tool_name)
    tool_layout.layout_main(self)

  def _show_error(self, error_message):
    ''' Show a pop-up dialog and print the error message.''' 
    # https://www.tutorialspoint.com/pyqt/pyqt_qmessagebox.htm
    print("TEMPY Error: " + error_message)
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("TEMPY Error")
    msg.setInformativeText(error_message)
    msg.setWindowTitle("Error")
    msg.setStandardButtons(QMessageBox.Cancel)
    msg.exec_()

  def _sccc_score(self):
    ''' Run the sccc score as a graphical function, 
    setting the colours of the chosen model.'''
    tool_scores.sccc_score(self)

  def _smoc_score(self):
    ''' Compute the smoc score but also plot
    the scores below the tool.'''
    tool_scores.smoc_score(self)

  def _nmi_score(self):
    ''' Run the nmi score, printing to the log'''
    tool_scores.nmi_score(self) 
  
  def onclick(self, event):
    ''' This is called when we click on the graph. We check to see which residue we are closest to by rounding down. '''
    # TODO - move this to a more appropriate place I think
    if event.x != None and event.y != None and event.xdata != None and event.ydata != None:
      #print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % ('double' if event.dblclick else 'single', event.button, event.x, event.y, event.xdata, event.ydata))
      res = int(math.floor(event.xdata))
       # Now find which line we are nearest
      sc = event.ydata 
      tc = math.fabs(float(self._current_smoc_scores[0][res][0]) - event.ydata)
      smodel = 0
      idx = 0

      for mm in self._current_smoc_scores:
        tb = math.fabs(float(mm[res][0]) - event.ydata) 
        if tb < tc:
          tc = tb
          smodel = idx
        idx += 1

      print("Clicked residue " + str(res) + " on model " + self._scored_models_smoc[smodel].__str__())
      # clear selection and reselect
      self.session.selection.clear()
      sres = self._scored_models_smoc[smodel].residues[res]
      
      for atom in sres.atoms:
        atom.selected = True
      # TODO - do we want bonds as well or are they already selected?

  def _select_rigid_file_sccc(self):
    from PyQt5.QtWidgets import QFileDialog
    filename = QFileDialog.getOpenFileName(None, 'OpenFile')
    print(filename)
    self._widget_rigid_file_sccc.setText(filename[0])

  def _select_rigid_file_smoc(self):
    from PyQt5.QtWidgets import QFileDialog
    filename = QFileDialog.getOpenFileName(None, 'OpenFile')
    print(filename)
    self._widget_rigid_file_smoc.setText(filename[0])


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
