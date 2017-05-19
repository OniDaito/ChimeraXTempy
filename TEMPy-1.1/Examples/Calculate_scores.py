from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
import os,sys
from TEMPy.class_arg import TempyParser
from traceback import print_exc

EXAMPLEDIR = 'Test_Files'
print 'use --help for help'
print '-m/-m1 [map] for input map, -m1,-m2 for two input maps'
print '-p/-p1 [pdb] for input pdb'
print '-r [resolution]; -r1,r2 for the two map resolutions'
print '-t [density threhold]; -t1,t2 for the two map thresholds'
print '\n\n###################################'

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
#EXAMPLE RUN
flag_example = False
if len(sys.argv) == 1:
  path_example=os.path.join(os.getcwd(),EXAMPLEDIR)
  if os.path.exists(path_example)==True:
    print "%s exists" %path_example
  #else: sys.exit('No input')
  #os.chdir(path_out)
  flag_example = True

#calculate map contour
def map_contour(m,t=-1.):
  mName = os.path.basename(m).split('.')[0]
  #print 'reading map'
  emmap=MapParser.readMRC(m)
  c1 = None
  if t != -1.0:
    print 'calculating contour'
    zeropeak,ave,sigma1 = emmap._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return mName,emmap,c1
#calculate model contour
def model_contour(p,res=4.0,emmap=False,t=-1.):
  pName,modelmap = blur_model(p,res,emmap)
  c1 = None
  if t != -1.0:
    print 'calculating contour'
    c1 = t*emmap.std()#0.0
  return pName,modelmap,c1
def blur_model(p,res=4.0,emmap=False):
  pName = os.path.basename(p).split('.')[0]
  print 'reading the model'
  structure_instance=PDBParser.read_PDB_file(pName,p,hetatm=False,water=False)
  print 'filtering the model'
  blurrer = StructureBlurrer()
  if res is None: sys.exit('Map resolution required..')
  #emmap = blurrer.gaussian_blur(structure_instance, res,densMap=emmap_1,normalise=True)
  modelmap = blurrer.gaussian_blur_real_space(structure_instance, res,densMap=emmap,normalise=True) 
  return pName,modelmap

#GET INPUT DATA
if flag_example:
  
  p = os.path.join(path_example,'1J6Z.pdb')
  m = os.path.join(path_example,'emd_5168_monomer.mrc')
  res = 6.6
  Name1 = os.path.basename(m).split('.')[0]
  Name2 = os.path.basename(p).split('.')[0]
  emmap1=MapParser.readMRC(m)
  structure_instance=PDBParser.read_PDB_file(Name2,p,hetatm=False,water=False)
  blurrer = StructureBlurrer()
  emmap2 = blurrer.gaussian_blur(structure_instance, res,densMap=emmap1)
  c1 = 9.7
  c2 = 1.0
elif all(x is None for x in [m,m1,m2]):
    # for 2 models
    if None in [p1,p2]:
        sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    Name1,emmap1,c1 = model_contour(p1,res=4.0,emmap=False,t=0.5)
    r1 = r2 = r = 4.0
    if c2 is None: Name2,emmap2,c2 = model_contour(p2,res=r,emmap=False,t=0.5)
    else: p2Name,emmap2 = blur_model(p2,res=r,emmap=False)
    flag_filt = False
    flag_scale = False
elif None in [m1,m2]:
    # for one map and model
    m = tp.args.inp_map
    print 'reading map'
    if c is None: Name1,emmap1,c1 = map_contour(m,t=1.5)
    else:
        Name1 = os.path.basename(m).split('.')[0]
        emmap1=MapParser.readMRC(m)
    if r1 is None and r is None: sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    elif r1 is None: r1 = r
    if all(x is None for x in [p,p1,p2]): sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    elif None in [p1,p2]: p = tp.args.pdb
    else: sys.exit('Input two maps or a map and model, map resolution(s) (required) and contours (optional)')
    r2 = 3.0
    #print 'reading model'
    Name2,emmap2,c2 = model_contour(p,res=r1,emmap=emmap1,t=0.5)
else: 
    # For 2 input maps
    if None in [r1,r2]: sys.exit('Input two maps, their resolutions(required) and contours(optional)')
    print 'reading map1'
    if c1 is None:
        Name1,emmap1,c1 = map_contour(m1,t=1.5)
    else:
        Name1 = os.path.basename(m1).split('.')[0]
        emmap1=MapParser.readMRC(m1)
    print 'reading map2' 
    if c2 is None:
        Name2,emmap2,c2 = map_contour(m2,t=1.5)
    else:
        Name2 = os.path.basename(m2).split('.')[0]
        emmap2=MapParser.readMRC(m2)
        



print 'Scoring...'
if not None in [Name1,Name2]:
  scores = {}
  sc = ScoringFunctions()
  #OVR
  try:
    ccc_mask,ovr = sc.CCC_map(emmap1,emmap2,c1,c2,3)
    print 'Percent overlap:', ovr
    if ovr < 0.0: ovr = 0.0
  except:
    print 'Exception for lccc and overlap score'
    print_exc()
    ovr = 0.0
  scores['overlap'] = ovr
  if ovr < 0.02:
    sys.exit("Maps do not overlap.")
  #SCCC
  print 'Local correlation score: ', ccc_mask
  if ccc_mask < -1.0 or ccc_mask > 1.0:
    ccc_mask = 0.0
  scores['local_correlation'] = ccc_mask
  #LMI
  try:
    mi_mask = sc.MI(emmap1,emmap2,c1,c2,3)
    print 'Local Mutual information score: ', mi_mask
    if mi_mask < 0.0: mi_mask = 0.0
  except:
    print 'Exception for MI score'
    print_exc()
    mi_mask = 0.0
  scores['mutual_information'] = mi_mask

  #NMI
  try:
    nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
    print 'Normalized Mutual information score:', nmi
    if nmi < 0.0: nmi = 0.0
  except:
    print 'Exception for NMI score'
    print_exc()
    nmi = 0.0
  #CD
  try:
    chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Minimum')
    if chm == 0.0 or chm is None:
      chm = sc._surface_distance_score(emmap1,emmap2,c1,c2,'Mean')
    print 'Surface distance score: ', chm
    if chm < 0.0: chm = 0.0
  except:
    print 'Exception for surface distance score'
    print_exc()
    chm = 0.0
  scores['surface_distance'] = chm
  
  try:
    #nv = sc.normal_vector_score(emmap_1,emmap_2,float(c1)-(emmap_1.std()*0.05),float(c1)+(emmap_1.std()*0.05),None)
    nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05),'Minimum')
    if nv == 0.0 or nv is None: nv = sc.normal_vector_score(emmap1,emmap2,float(c1),float(c1)+(emmap1.std()*0.05))
    print 'Normal vector score: ', nv
    if nv < 0.0: 
      nv = 0.0
  except:
    print 'Exception for NV score'
    print_exc()
    nv = 0.0
  scores['normal_vector'] = nv
