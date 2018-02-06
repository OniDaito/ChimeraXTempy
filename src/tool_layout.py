# vim: set expandtab shiftwidth=2 softtabstop=2:

# tool_layout.py
# An additional file to perform the QT layout a little more neatly

from PyQt5.QtWidgets import QLineEdit, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QPushButton, QWidget, QTabWidget, QSizePolicy, QCheckBox, QListWidget

from PyQt5 import QtCore

from chimerax.core.ui.gui import MainToolWindow

def _adv_clicked(tool):
  ''' A closure for our advanced options checkbox.'''
  def _on_clicked(state):
    if state == QtCore.Qt.Checked:
      tool._widget_c1_nmi.setEnabled(True)
      tool._widget_c2_nmi.setEnabled(True)
    else:
      tool._widget_c1_nmi.setEnabled(False)
      tool._widget_c2_nmi.setEnabled(False)  
  return _on_clicked

def _adv_dif_clicked(tool):
  ''' A closure for our advanced options checkbox.'''
  def _on_clicked(state):
    if state == QtCore.Qt.Checked:
      tool._widget_c1_dif.setEnabled(True)
      tool._widget_c2_dif.setEnabled(True)
    else:
      tool._widget_c1_dif.setEnabled(False)
      tool._widget_c2_dif.setEnabled(False)  
  return _on_clicked


def layout_nmi(tool, tab):
  ''' Layout the nmi tab '''
  layout = QVBoxLayout()
  tool.nmi_layout = layout
  layout.setContentsMargins(0, 5, 0, 5)
  tab.setLayout(layout)
  layout.setSpacing(5)
  layout.setAlignment(QtCore.Qt.AlignTop)
 
  params_layout = QHBoxLayout()
  layout.addLayout(params_layout)
  params_layout2 = QHBoxLayout()

 
  label_rez1 = QLabel("Resolution.1")
  label_rez2 = QLabel("Resolution.2")
  label_c1 = QLabel("contour.1")
  label_c2 = QLabel("contour.2")
  
  label_rez1.setFixedSize(120,30)
  label_rez2.setFixedSize(120,30)
  label_c1.setFixedSize(120,30)
  label_c2.setFixedSize(120,30)
  
  tool._widget_rez1_nmi = QLineEdit()
  tool._widget_rez2_nmi = QLineEdit()
  tool._widget_c1_nmi = QLineEdit()
  tool._widget_c2_nmi = QLineEdit()
  
  adv_box = QCheckBox("Advanced Options")
  ac = _adv_clicked(tool)
  adv_box.stateChanged.connect(ac)
  tool._widget_c1_nmi.setEnabled(False)
  tool._widget_c2_nmi.setEnabled(False)

  tool._widget_rez1_nmi.setFixedSize(40,30)
  tool._widget_rez2_nmi.setFixedSize(40,30)
  tool._widget_c1_nmi.setFixedSize(40,30)
  tool._widget_c2_nmi.setFixedSize(40,30)

  params_layout.addWidget(label_rez1)
  params_layout.addWidget(tool._widget_rez1_nmi)
  params_layout.addWidget(label_rez2)
  params_layout.addWidget(tool._widget_rez2_nmi)
 
  layout.addWidget(adv_box)

  layout.addLayout(params_layout2)
  params_layout2.addWidget(label_c1)
  params_layout2.addWidget(tool._widget_c1_nmi)
  params_layout2.addWidget(label_c2)
  params_layout2.addWidget(tool._widget_c2_nmi)

  tool._widget_rez1_nmi.setText("10.0")
  tool._widget_rez2_nmi.setText("10.0")
  tool._widget_c1_nmi.setText("")
  tool._widget_c2_nmi.setText("")
 
  button_nmi = QPushButton("NMI")
  button_nmi.clicked.connect(tool._nmi_score)
  tool.nmi_layout.addWidget(button_nmi)

def layout_dif(tool, tab):
  ''' Layout the difmap tab '''
  layout = QVBoxLayout()
  tool.dif_layout = layout
  layout.setContentsMargins(0, 5, 0, 5)
  tab.setLayout(layout)
  layout.setSpacing(5)
  layout.setAlignment(QtCore.Qt.AlignTop)
 
  params_layout = QHBoxLayout()
  layout.addLayout(params_layout)
  params_layout2 = QHBoxLayout()

  label_rez1 = QLabel("Resolution.1")
  label_rez2 = QLabel("Resolution.2")
  label_c1 = QLabel("contour.1")
  label_c2 = QLabel("contour.2")
  
  label_rez1.setFixedSize(120,30)
  label_rez2.setFixedSize(120,30)
  label_c1.setFixedSize(120,30)
  label_c2.setFixedSize(120,30)
  
  tool._widget_rez1_dif = QLineEdit()
  tool._widget_rez2_dif = QLineEdit()
  tool._widget_c1_dif = QLineEdit()
  tool._widget_c2_dif = QLineEdit()
  
  adv_box = QCheckBox("Advanced Options")
  ac = _adv_dif_clicked(tool)
  adv_box.stateChanged.connect(ac)
  tool._widget_c1_dif.setEnabled(False)
  tool._widget_c2_dif.setEnabled(False)

  tool._widget_rez1_dif.setFixedSize(40,30)
  tool._widget_rez2_dif.setFixedSize(40,30)
  tool._widget_c1_dif.setFixedSize(40,30)
  tool._widget_c2_dif.setFixedSize(40,30)

  params_layout.addWidget(label_rez1)
  params_layout.addWidget(tool._widget_rez1_dif)
  params_layout.addWidget(label_rez2)
  params_layout.addWidget(tool._widget_rez2_dif)
 
  layout.addWidget(adv_box)

  layout.addLayout(params_layout2)
  params_layout2.addWidget(label_c1)
  params_layout2.addWidget(tool._widget_c1_dif)
  params_layout2.addWidget(label_c2)
  params_layout2.addWidget(tool._widget_c2_dif)

  tool._widget_rez1_dif.setText("10.0")
  tool._widget_rez2_dif.setText("10.0")
  tool._widget_c1_dif.setText("")
  tool._widget_c2_dif.setText("")
 
  button_dif = QPushButton("Generate DifMap")
  button_dif.clicked.connect(tool._gen_difmap)
  tool.dif_layout.addWidget(button_dif)

def _sccc_scores_clicked(toptool):
  ''' A little closure over the tool in question.'''
  def clicked(item):
    idx = toptool._sccc_scores_widget.currentRow()
    RB, score = toptool._current_sccc_scores[idx]
    atomic_model = toptool._scored_models_sccc
    
    toptool.session.selection.clear()
    for atom in RB.get_atom_list():
      sres = atomic_model.residues[atom.res_no]
      # Duplicate work here but who cares
      for atom in sres.atoms:
        atom.selected = True
   
  return clicked

def add_sccc_scores(tool, scores):
  ''' When we have scores, clear existing and add new ones.'''
  tool._sccc_scores_widget.clear()
  for RB,score in scores:
    tool._sccc_scores_widget.addItem(str(RB) + ": " + str(score))

def layout_sccc(tool, tab):
  ''' Layout the sccc tab '''
  layout = QVBoxLayout()
  tool.sccc_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)
  layout.setSpacing(5)
  layout.setAlignment(QtCore.Qt.AlignTop)
  
  params_layout = QHBoxLayout()
  params_layout.setSpacing(0)
  layout.addLayout(params_layout)
  params_layout.setContentsMargins(0, 5, 0, 5)
  
  label_rez = QLabel("Resolution")
  label_sigma = QLabel("Sigma")

  label_rez.setFixedSize(90,30)
  label_sigma.setFixedSize(60,30)
  
  tool._widget_sigma_sccc = QLineEdit()
  tool._widget_rez_sccc = QLineEdit()
  tool._widget_window_sccc = QLineEdit()

  tool._widget_rez_sccc.setFixedSize(40,30)
  tool._widget_sigma_sccc.setFixedSize(50,30)

  params_layout.addWidget(label_rez)
  params_layout.addWidget(tool._widget_rez_sccc)
  params_layout.addWidget(label_sigma)
  params_layout.addWidget(tool._widget_sigma_sccc)

  tool._widget_rez_sccc.setText("10.0")
  tool._widget_sigma_sccc.setText("0.187")
  tool._widget_window_sccc.setText("9")

  # Rigid file
  rigid_layout = QHBoxLayout()
  layout.addLayout(rigid_layout)
  button_file = QPushButton("Rigid file")
  tool._widget_rigid_file_sccc = QLineEdit()
  # TODO - remove this eventually
  tool._widget_rigid_file_sccc.setText('/home/oni/Projects/ChimeraXTempy/test/rigid_RF.txt')

  button_file.clicked.connect(tool._select_rigid_file_sccc)
  rigid_layout.addWidget(button_file)
  rigid_layout.addWidget(tool._widget_rigid_file_sccc)

  button_sccc = QPushButton("SCCC")
  button_sccc.clicked.connect(tool._sccc_score)
  tool.sccc_layout.addWidget(button_sccc)

  list_scores = QListWidget()
  tool.sccc_layout.addWidget(list_scores)
  tool._sccc_scores_widget = list_scores
  list_scores.show()
  list_scores.clicked.connect(_sccc_scores_clicked(tool))

def layout_smoc(tool, tab):
  ''' Layout the smoc tab '''
  layout = QVBoxLayout()
  tool.smoc_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)
  layout.setSpacing(5)
  layout.setAlignment(QtCore.Qt.AlignTop)
 
  params_layout = QHBoxLayout()
  params_layout.setContentsMargins(0, 5, 0, 5)
  layout.addLayout(params_layout)
  
  label_rez = QLabel("Resolution.")
  label_sigma = QLabel("Sigma")
  label_window = QLabel("Window")

  label_rez.setFixedSize(90,30)
  label_sigma.setFixedSize(60,30)
  label_window.setFixedSize(70,30)

  tool._widget_sigma_smoc = QLineEdit()
  tool._widget_rez_smoc = QLineEdit()
  tool._widget_window_smoc = QLineEdit()

  tool._widget_rez_smoc.setFixedSize(40,30)
  tool._widget_sigma_smoc.setFixedSize(50,30)
  tool._widget_window_smoc.setFixedSize(40,30)

  params_layout.addWidget(label_rez)
  params_layout.addWidget(tool._widget_rez_smoc)
  params_layout.addWidget(label_sigma)
  params_layout.addWidget(tool._widget_sigma_smoc)
  params_layout.addWidget(label_window)
  params_layout.addWidget(tool._widget_window_smoc)

  tool._widget_rez_smoc.setText("10.0")
  tool._widget_sigma_smoc.setText("0.187")
  tool._widget_window_smoc.setText("9")

  # Rigid file
  rigid_layout = QHBoxLayout()
  layout.addLayout(rigid_layout)
  button_file = QPushButton("Rigid file")
  tool._widget_rigid_file_smoc = QLineEdit()
  # TODO - remove this eventually
  tool._widget_rigid_file_smoc.setText('/home/oni/Projects/ChimeraXTempy/test/rigid_RF.txt')
  button_file.clicked.connect(tool._select_rigid_file_smoc)
  rigid_layout.addWidget(button_file)
  rigid_layout.addWidget(tool._widget_rigid_file_smoc)

  button_smoc = QPushButton("SMOC")
  button_smoc.clicked.connect(tool._smoc_score)
  tool.smoc_layout.addWidget(button_smoc)

def layout_ccc(tool, tab):
  ''' Layout the ccc tab '''
  layout = QVBoxLayout()
  tool.ccc_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)
  layout.setSpacing(5)
  layout.setAlignment(QtCore.Qt.AlignTop)

  params_layout = QHBoxLayout()
  params_layout.setContentsMargins(0, 5, 0, 5)
  layout.addLayout(params_layout)
  
  label_rez = QLabel("Resolution.")
  label_rez.setFixedSize(90,30)

  tool._widget_rez_ccc = QLineEdit()
  tool._widget_rez_ccc.setFixedSize(40,30)
  
  params_layout.addWidget(label_rez)
  params_layout.addWidget(tool._widget_rez_ccc)

  tool._widget_rez_smoc.setText("10.0")

  button_ccc = QPushButton("CCC")
  button_ccc.clicked.connect(tool._ccc_score)
  tool.ccc_layout.addWidget(button_ccc)

def layout_main(tool):
  # Standard template stuff 
  tool.display_name = "Tempy"
  tool.tool_window = MainToolWindow(tool)
  tool.tool_window.manage(placement="side")
  parent = tool.tool_window.ui_area
  policy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
  parent.setSizePolicy(policy)
  layout = QVBoxLayout()
  tool.top_layout = layout

  # Create some tabs
  tab_widget = QTabWidget(parent)
  tab_widget.setMinimumHeight(210) 
  tab_widget.setMinimumWidth(400)
  tab_widget.setSizePolicy(policy)

  tab_smoc = QWidget()
  tab_smoc.setMinimumHeight(210) 
  tab_smoc.setMinimumWidth(400)
  tab_smoc.setSizePolicy(policy)
  tab_widget.addTab(tab_smoc, "SMOC")

  tab_sccc = QWidget()
  tab_sccc.setMinimumHeight(210) 
  tab_sccc.setMinimumWidth(400)
  tab_sccc.setSizePolicy(policy)
  tab_widget.addTab(tab_sccc, "SCCC")

  tab_ccc = QWidget()
  tab_ccc.setMinimumHeight(210) 
  tab_ccc.setMinimumWidth(400)
  tab_ccc.setSizePolicy(policy)
  tab_widget.addTab(tab_ccc, "CCC")

  tab_nmi = QWidget()
  tab_nmi.setMinimumHeight(210) 
  tab_nmi.setMinimumWidth(400)
  tab_nmi.setSizePolicy(policy)
  tab_widget.addTab(tab_nmi, "NMI")
  
  tab_dif = QWidget()
  tab_dif.setMinimumHeight(210) 
  tab_dif.setMinimumWidth(400)
  tab_dif.setSizePolicy(policy)
  tab_widget.addTab(tab_dif, "DIFMAP")

  # Now layout each tab
  layout_nmi(tool, tab_nmi)
  layout_sccc(tool, tab_sccc)
  layout_smoc(tool, tab_smoc)
  layout_ccc(tool, tab_ccc)
  layout_dif(tool, tab_dif)

  parent.setLayout(layout)
  layout.addWidget(tab_widget)

  # Figures for SMOC
  tool._figure = None
  tool._canvas = None
