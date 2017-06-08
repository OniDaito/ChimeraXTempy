# vim: set expandtab shiftwidth=2 softtabstop=2:

# sccc
#   "score" - compute the sccc score using tempy
#

from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser
from TEMPy.ProtRep_Biopy import BioPy_Structure,BioPyAtom

from .util import chimera_to_tempy_atom, chimera_to_tempy_map

def score(session, atomic_model, map_model, rigid_filename):
    print("Calculating SCCC Score")
    rez = 10.0

   
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
    listRB = RBParser.read_FlexEM_RIBFIND_files(rigid_filename, bio_atom_structure)
    
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

