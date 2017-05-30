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

from chimerax.core.tools import ToolInstance

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
    from PyQt5.QtWidgets import QGridLayout, QLabel, QComboBox, QPushButton
    from chimerax.core.ui.widgets import HtmlView
    layout = QGridLayout()
    layout.setContentsMargins(0, 0, 0, 0)
    label = QLabel("SCCC Score:")
    layout.addWidget(label, 0, 0)
    button = QPushButton("calculate")
    button.clicked.connect(self._sccc_score)
    layout.addWidget(button, 0, 2)
    self.results_view = HtmlView(parent, size_hint=(575, 300),
                                 interceptor=self._html_view,
                                 schemes=[])
    layout.addWidget(self.results_view, 1, 0, 1, 3)
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
    tp = TempyParser()
    # OPTIONS
    m = "filemap"#tp.args.inp_map - current mrc file apparently
    r = 10.0 #tp.args.res
    c = "dunno" #tp.args.thr
    p = "thepdb" #tp.args.pdb - current pdb file
    
    # the sigma factor determines the width of the Gaussian distribution used to describe each atom
    sim_sigma_coeff = 0.187
    if not tp.args.sigfac is None:
        sim_sigma_coeff = tp.args.sigfac
    
    return
    #GET INPUT DATA
    if flag_example:
        m = os.path.join(path_example,'1akeA_10A.mrc')
        p = os.path.join(path_example,'1ake_mdl1.pdb')
        r = 10.0
        rb_file = os.path.join(path_example,'1ake_mdl1_rigid.txt')
    elif None in [m1,m2]:
        # for one map and model
        m = tp.args.inp_map
        print('reading map')
        Name1 = os.path.basename(m).split('.')[0]
        emmap1=MapParser.readMRC(m)
        if r1 is None and r is None: sys.exit('Input a map and model, map resolution (required)')
        elif r1 is None: r1 = r
        if all(x is None for x in [p,p1,p2]): sys.exit('Input a map and model, map resolution (required)')
        elif None in [p1,p2]: p = tp.args.pdb
        else: sys.exit('Input a map and model, map resolution (required)')
        rb_file = tp.args.rigidfile
        if rb_file is None: sys.exit('Rigid body file missing')

    # make class instances for density simulation (blurring), scoring and plot scores
    blurrer = StructureBlurrer()
    scorer = ScoringFunctions()

    # read map file
    emmap=MapParser.readMRC(m)
    # read PDB file
    structure_instance=PDBParser.read_PDB_file('pdbfile',p,hetatm=False,water=False)
   
    SCCC_list_structure_instance=[]
    # read rigid body file and generate structure instances for each segment
    listRB=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance)
    # score each rigid body segment
    listsc_sccc = []
    print('calculating scores')
    for RB in listRB:
      # sccc score
      score_SCCC=scorer.SCCC(emmap,r,sim_sigma_coeff,structure_instance,RB,c_mode=False)
      SCCC_list_structure_instance.append(score_SCCC)
      print ('>>', score_SCCC)
      listsc_sccc.append(score_SCCC)
      
    listRB=RBParser.RBfileToRBlist(rb_file)
    if len(listRB) == len(listsc_sccc):
      #include sccc scores as b-factor records
      for x in structure_instance.atomList:
        cur_chain = x.chain
        cur_res = x.get_res_no()
        ct = 0
        flage = 0
        if cur_chain in ['',' ']: cur_chain = '-'
        for rb in listRB:
          for rb1 in rb:
            if len(rb1) == 2:
                try: st = int(rb1[0].split(':')[0])
                except TypeError: st = rb1[0]
                try: en = int(rb1[1].split(':')[0])
                except TypeError: en = rb1[1]
                if ':' in rb1[0]:
                  ch_rb = rb1[0].split(':')[1]
                else: ch_rb = '-'
                #TODO check for insertion codes
                for i in range(st,en+1):
                  if int(cur_res) == i and ch_rb == cur_chain: 
                    flage = 1
                    break
                if flage == 1: break
          if flage == 1: break  
          ct += 1
        if flage == 1:
          sc = listsc_sccc[ct]    
          try: x.temp_fac = sc
          except :
            print('Residue missing: ',cur_res, ch, out_iter_pdb)
            x.temp_fac = 0.0
        else:
          x.temp_fac = 0.0
      pName = os.path.basename(os.path.abspath(p)).split('.')[0]
      structure_instance.write_to_PDB(os.path.join(os.path.dirname(p),pName+"_sc.pdb"))  

    # generate chimera attribute file for coloring segments based on sccc score
    # Plot.PrintOutChimeraAttributeFileSCCC_Score(p,SCCC_list_structure_instance,listRB)
    # TODO - do this directly

    #if os.path.isfile(os.path.abspath(p)):
    #  pName = os.path.basename(os.path.abspath(p)).split('.')[0]
    #  scf = open(os.path.join(os.path.dirname(os.path.abspath(p)),'sccc_'+pName),'w')
    #  for sc in listsc_sccc: scf.write(str(sc)+"\n")
    #  scf.close()
      


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
