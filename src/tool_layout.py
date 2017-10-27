# vim: set expandtab shiftwidth=2 softtabstop=2:

# tool_layout.py
# An additional file to perform the QT layout a little more neatly

from PyQt5.QtWidgets import QLineEdit, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QPushButton, QWidget, QTabWidget, QSizePolicy

from chimerax.core.ui.gui import MainToolWindow

def layout_nmi(tool, tab):
  ''' Layout the nmi tab '''
  layout = QVBoxLayout()
  tool.nmi_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)

  button_nmi = QPushButton("NMI")
  button_nmi.clicked.connect(tool._nmi_score)
  tool.nmi_layout.addWidget(button_nmi)
  
  label_rez = QLabel("Res.")
  label_rez2 = QLabel("Res.2")
  label_c1 = QLabel("c.1")
  label_c2 = QLabel("c.2")
  
  label_rez.setFixedSize(40,30)
  label_rez2.setFixedSize(40,30)
  label_c1.setFixedSize(40,30)
  label_c2.setFixedSize(40,30)
  
  tool._widget_rez = QLineEdit()
  tool._widget_rez2 = QLineEdit()
  tool._widget_c1 = QLineEdit()
  tool._widget_c2 = QLineEdit()

  tool._widget_rez.setFixedSize(40,30)
  tool._widget_rez2.setFixedSize(40,30)
  tool._widget_c1.setFixedSize(40,30)
  tool._widget_c2.setFixedSize(40,30)

  layout.addWidget(label_rez)
  layout.addWidget(tool._widget_rez)
  layout.addWidget(label_rez2)
  layout.addWidget(tool._widget_rez2)
  layout.addWidget(label_c1)
  layout.addWidget(tool._widget_c1)
  layout.addWidget(label_c2)
  layout.addWidget(tool._widget_c2)

  tool._widget_rez.setText("10.0")
  tool._widget_rez2.setText("10.0")
  tool._widget_c1.setText("10.0")
  tool._widget_c2.setText("10.0")

def layout_sccc(tool, tab):
  ''' Layout the sccc tab '''
  layout = QVBoxLayout()
  tool.sccc_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)

  button_sccc = QPushButton("SCCC")
  button_sccc.clicked.connect(tool._sccc_score)
  tool.sccc_layout.addWidget(button_sccc)

  label_rez = QLabel("Res.")
  label_sigma = QLabel("Sigma")

  label_rez.setFixedSize(40,30)
  label_sigma.setFixedSize(60,30)
  
  tool._widget_sigma_sccc = QLineEdit()
  tool._widget_rez_sccc = QLineEdit()
  tool._widget_window_sccc = QLineEdit()

  tool._widget_rez_sccc.setFixedSize(40,30)
  tool._widget_sigma_sccc.setFixedSize(50,30)

  layout.addWidget(label_rez)
  layout.addWidget(tool._widget_rez_sccc)
  layout.addWidget(label_sigma)
  layout.addWidget(tool._widget_sigma_sccc)

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

def layout_smoc(tool, tab):
  ''' Layout the smoc tab '''
  layout = QVBoxLayout()
  tool.smoc_layout = layout
  layout.setContentsMargins(0, 0, 0, 0)
  tab.setLayout(layout)

  button_smoc = QPushButton("SMOC")
  button_smoc.clicked.connect(tool._smoc_score)
  tool.smoc_layout.addWidget(button_smoc)

  label_rez = QLabel("Res.")
  label_sigma = QLabel("Sigma")
  label_window = QLabel("Window")

  label_rez.setFixedSize(40,30)
  label_sigma.setFixedSize(60,30)
  label_window.setFixedSize(70,30)

  tool._widget_sigma_smoc = QLineEdit()
  tool._widget_rez_smoc = QLineEdit()
  tool._widget_window_smoc = QLineEdit()

  tool._widget_rez_smoc.setFixedSize(40,30)
  tool._widget_sigma_smoc.setFixedSize(50,30)
  tool._widget_window_smoc.setFixedSize(40,30)

  layout.addWidget(label_rez)
  layout.addWidget(tool._widget_rez_smoc)
  layout.addWidget(label_sigma)
  layout.addWidget(tool._widget_sigma_smoc)
  layout.addWidget(label_window)
  layout.addWidget(tool._widget_window_smoc)

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


def layout_main(tool):
  # Standard template stuff 
  tool.display_name = "Tempy"
  tool.tool_window = MainToolWindow(tool)
  tool.tool_window.manage(placement="side")
  parent = tool.tool_window.ui_area
  policy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
  tool.setMinimumHeight(300)
  
  layout = QVBoxLayout()

  # Create some tabs
  tab_widget = QTabWidget(parent)
  tab_widget.setMinimumHeight(300) 
  tab_widget.setMinimumWidth(400)
  tab_widget.setSizePolicy(policy)

  tab_smoc = QWidget()
  tab_smoc.setMinimumHeight(300) 
  tab_smoc.setMinimumWidth(400)
  tab_smoc.setSizePolicy(policy)
  tab_widget.addTab(tab_smoc, "SMOC")

  tab_sccc = QWidget()
  tab_sccc.setMinimumHeight(300) 
  tab_sccc.setMinimumWidth(400)
  tab_sccc.setSizePolicy(policy)
  tab_widget.addTab(tab_sccc, "SCCC")

  tab_nmi = QWidget()
  tab_nmi.setMinimumHeight(300) 
  tab_nmi.setMinimumWidth(400)
  tab_nmi.setSizePolicy(policy)
  tab_widget.addTab(tab_nmi, "NMI")

  # Now layout each tab
  layout_nmi(tool, tab_nmi)
  layout_sccc(tool, tab_sccc)
  layout_smoc(tool, tab_smoc)

  parent.setLayout(layout)
  # Figures for SMOC
  tool._figure = None
  tool._canvas = None
