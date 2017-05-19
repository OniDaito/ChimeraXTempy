import sys,os
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
from TEMPy.RigidBodyParser import RBParser
from TEMPy.ShowPlot import Plot
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.class_arg import TempyParser

#from datetime import datetime
EXAMPLEDIR = 'Test_Files'
'''
print 'use --help for help'
print '-m/-m1 [map] for input map'
print '-p/-p1 [pdb] for input pdb'
print '-r [resolution]'
print '-rf [rigid body file]'
print '-t [density threhold], map density threshold (default: 1.5*sigma)'
print '--sf [sigma factor], determines the width of the Gaussian used to describe each atom (default 0.187; other values 0.225,0.356,0.425)'

print '\n\n###################################'
'''
tp = TempyParser()
tp.generate_args()
#COMMAND LINE OPTIONS
m1 = tp.args.inp_map1
m2 = tp.args.inp_map2
m = tp.args.inp_map
r1 = tp.args.res1
r2 = tp.args.res2
r = tp.args.res
c1 = tp.args.thr1
c2 = tp.args.thr2
c = tp.args.thr
p = tp.args.pdb
p1 = tp.args.pdb1
p2 = tp.args.pdb2
#apix = tp.args.apix
# the sigma factor determines the width of the Gaussian distribution used to describe each atom
sim_sigma_coeff = 0.187
if not tp.args.sigfac is None:
    sim_sigma_coeff = tp.args.sigfac
#EXAMPLE RUN
flag_example = False
if len(sys.argv) == 1:
  path_example=os.path.join(os.getcwd(),EXAMPLEDIR)
  if os.path.exists(path_example)==True:
    print "%s exists" %path_example
  #else: sys.exit('No input')
  #os.chdir(path_out)
  flag_example = True
  
#GET INPUT DATA
if flag_example:
    m = os.path.join(path_example,'1akeA_10A.mrc')
    p = os.path.join(path_example,'1ake_mdl1.pdb')
    r = 10.0
    rb_file = os.path.join(path_example,'1ake_mdl1_rigid.txt')
elif None in [m1,m2]:
    # for one map and model
    m = tp.args.inp_map
    print 'reading map'
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
Plot=Plot()

# read map file
emmap=MapParser.readMRC(m)
# read PDB file
structure_instance=PDBParser.read_PDB_file('pdbfile',p,hetatm=False,water=False)
# generate atom density and blur to required resolution
#sim_map = blurrer.gaussian_blur(structure_instance, r,densMap=emmap,sigma_coeff=sim_sigma_coeff,normalise=True)

#sim_map = blurrer.gaussian_blur_real_space(structure_instance, r,densMap=emmap,sigma_coeff=sim_sigma_coeff,normalise=True)
SCCC_list_structure_instance=[]
# read rigid body file and generate structure instances for each segment
listRB=RBParser.read_FlexEM_RIBFIND_files(rb_file,structure_instance)
# score each rigid body segment
listsc_sccc = []
print 'calculating scores'
for RB in listRB:
  # sccc score
  score_SCCC=scorer.SCCC(emmap,r,sim_sigma_coeff,structure_instance,RB)
  SCCC_list_structure_instance.append(score_SCCC)
  print '>>', score_SCCC
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
      '''
      if len(rb) == 2:
        try: st = int(rb[0].split(':')[0])
        except TypeError: st = rb[0]
        try: en = int(rb[1].split(':')[0])
        except TypeError: en = rb[1]
        #TODO check for insertion codes
        for i in range(st,en):
          if int(cur_res) == i: 
            flage = 1
            break
      '''
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
      except KeyError, IndexError:
        print 'Residue missing: ',cur_res, ch, out_iter_pdb
        x.temp_fac = 0.0
    else:
      x.temp_fac = 0.0
  pName = os.path.basename(os.path.abspath(p)).split('.')[0]
  structure_instance.write_to_PDB(os.path.join(os.path.dirname(p),pName+"_sc.pdb"))  

# generate chimera attribute file for coloring segments based on sccc score
Plot.PrintOutChimeraAttributeFileSCCC_Score(p,SCCC_list_structure_instance,listRB)

if os.path.isfile(os.path.abspath(p)):
  pName = os.path.basename(os.path.abspath(p)).split('.')[0]
  scf = open(os.path.join(os.path.dirname(os.path.abspath(p)),'sccc_'+pName),'w')
  for sc in listsc_sccc: scf.write(str(sc)+"\n")
  scf.close()
  

