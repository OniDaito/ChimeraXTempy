# vim: set expandtab shiftwidth=2 softtabstop=2:

# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
# Seems like Tool exists inside the bundle and needs overridding


from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom

from chimerax.core.tools import ToolInstance
from chimerax.core.models import Models

from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

from .util import chimera_to_tempy_atom, chimera_to_tempy_map

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
    self._widget_rigid_file.setText('/home/oni/Projects/tempy/data/tempy/rigid.txt')

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
    print("Calculating SMOC Score")
    rez = 10.0 #tp.args.res
  
    atomic_model = None
    map_model = None

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
    
    # the sigma factor determines the width of the Gaussian distribution used to describe each atom
    # TODO - this needs to be a slider or textbox
    sim_sigma_coeff = 0.187

 
    # TODO - replace this with the current loaded map
    #emmap = MapParser.readMRC(map_model_filename)

    # make class instances for density simulation (blurring), scoring and plot scores
    blurrer = StructureBlurrer()
    scorer = ScoringFunctions()

    # read map file
    #emmap=MapParser.readMRC(m)
    
    # read PDB file
    # TODO - we need to change this to current selected PDB Model
    #structure_instance=PDBParser.read_PDB_file('pdbfile',p,hetatm=False,water=False)
  
    atomlist = []

    for atom in atomic_model.atoms:
      atomlist.append(chimera_to_tempy_atom(atom, len(atomlist)))

    bio_atom_structure = BioPy_Structure(atomlist)
    bio_map_structure = chimera_to_tempy_map(map_model)
    SCCC_list_structure_instance=[]
    # read rigid body file and generate structure instances for each segment
    listRB = RBParser.read_FlexEM_RIBFIND_files(rb_file, bio_atom_structure)
    
    # score each rigid body segment
    listsc_sccc = []
    print('calculating scores')

    for RB in listRB:
      # sccc score
      score_SCCC=scorer.SCCC(bio_map_structure, rez, sim_sigma_coeff, bio_atom_structure, RB, c_mode=False)
      SCCC_list_structure_instance.append(score_SCCC)
      print ('>>', score_SCCC)
      listsc_sccc.append(score_SCCC)
    
    # generate chimera attribute file for coloring segments based on sccc score
    # Plot.PrintOutChimeraAttributeFileSCCC_Score(p,SCCC_list_structure_instance,listRB)
    # TODO - do this directly

    #if os.path.isfile(os.path.abspath(p)):
    #  pName = os.path.basename(os.path.abspath(p)).split('.')[0]
    #  scf = open(os.path.join(os.path.dirname(os.path.abspath(p)),'sccc_'+pName),'w')
    #  for sc in listsc_sccc: scf.write(str(sc)+"\n")
    #  scf.close()
  
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
