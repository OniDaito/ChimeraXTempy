# vim: set expandtab shiftwidth=2 softtabstop=2:

# nmi
#   "score" - compute the nmi score using tempy
#

from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.class_arg import TempyParser

from traceback import print_exc
import os,sys

from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

from .util import chimera_to_tempy_model, chimera_to_tempy_map

#calculate map contour
def map_contour(m,t=-1.):
  c1 = None
  if t != -1.0:
    zeropeak,ave,sigma1 = m._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return c1

#calculate model contour
def model_contour(p,res=4.0,emmap=False,t=-1.):
  blurrer = StructureBlurrer()
  modelmap = blurrer.gaussian_blur_real_space(p, res,densMap=emmap,normalise=True) 
  c1 = None
  if t != -1.0:
    c1 = t*emmap.std()#0.0
  return modelmap,c1

def lpfilter(emmap,r):
  cutoff = emmap.apix/float(r)
  mapfilt = emmap._tanh_lowpass(cutoff)
  return mapfilt

def match_grid(emmap1,emmap2,c1,c2):
  # DETERMINE A COMMON ALIGNMENT BOX : fill minvalue for extra voxel pads
  spacing = emmap2.apix
  if emmap2.apix < emmap1.apix: spacing = emmap1.apix
  grid_shape, new_ori = emmap1._alignment_box(emmap2,spacing)
  # INTERPOLATE TO NEW GRID
  try: emmap_1 = emmap1._interpolate_to_grid1(grid_shape,spacing,new_ori)
  except: emmap_1 = emmap1._interpolate_to_grid(grid_shape,spacing,new_ori)
  try: c1 = emmap_1._find_level(np.sum(emmap1.fullMap>c1)*(emmap1.apix**3))
  except: pass
  del emmap1.fullMap
  del emmap1
  try: emmap_2 = emmap2._interpolate_to_grid1(grid_shape,spacing,new_ori)
  except: emmap_2 = emmap2._interpolate_to_grid(grid_shape,spacing,new_ori)
  try: c2 = emmap_2._find_level(np.sum(emmap2.fullMap>c2)*(emmap2.apix**3))
  except: pass
  del emmap2.fullMap
  del emmap2
  return emmap_1, emmap_2

def score_cmd(session, comparators, compared, rez_comparators, rez_compared, contours_comparators, contour_compared):
  sc = ScoringFunctions()
  blurrer = StructureBlurrer()

   # Loop through these to be compared
  idx = 0
  scores = []
  
  for comparator in comparators:
    emmap1 = None
    emmap2 = None
    
    if type(comparator) is AtomicStructure:
      if type(compared) is AtomicStructure:
        # Both models
        if None in ([rez_compared] + rez_comparators): 
          print("Please provide the resolution for all models")
          return
       
        bms1 = chimera_to_tempy_model(compared)
        bms2 = chimera_to_tempy_model(comparator)
        emmap1 = model_contour( bms1, rez_compared, emmap=False,t=0.5)
        
        if contours_comparators[idx] is None: 
          emmap2 = model_contour(bms2, rez_comparators[idx],emmap=False,t=0.5)
        else:
          emmap2 = blur_model(bms2, rez_comparators[idx], emmap=False)
     
      else:
        # 0 - map, 1 - model
        if rez_comparators[idx] == None: 
          print("Please provide the resolution for the model.")
          return

        emmap1 = chimera_to_tempy_map(compared)
        bms = chimera_to_tempy_model(comparator)
        emmap2 = blurrer.gaussian_blur(bms, rez_compared, densMap=emmap1)
        
    else:
      if type(compared) is AtomicStructure:
        # 0 - model, 1 - map
        if rez_compared == None: 
          print("Please provide the resolution for the model.")
          return

        emmap2 = chimera_to_tempy_map(comparator)
        bms = chimera_to_tempy_model(compared)
        emmap1 = blurrer.gaussian_blur(bms, rez_compared, densMap=emmap2)

      else:
        # 0 - map, 1 - map
        emmap1 = chimera_to_tempy_map(compared)
        emmap2 = chimera_to_tempy_map(comparator)
     
    c1 = contour_compared
    # Contouring
    if c1 == None:
      c1 = map_contour(emmap1,t=1.5)

    c2 = contours_comparators[idx]
    # This kinda makes no sense and could be tricky
    if c2 == None:
      c2 = map_contour(emmap2,t=1.5)

    # Some kind of fix if the maps don't match?
    # Resize, resample or blur of somekind
    if not sc.mapComparison(emmap1,emmap2):
      emmap1._crop_box(c1,0.5)
      emmap2._crop_box(c2,0.5)
      
      if rez_compared > 1.25*rez_comparators[idx]: 
        emmap_2 = lpfilter(emmap2,rez_compared)
        emmap1, emmap2 = match_grid(emmap1,emmap_2,c1,c2)
      elif rez_comparators[idx] > 1.25*rez_compared:
        emmap_1 = lpfilter(emmap1,rez_comparators[idx])
        emmap1, emmap2 = match_grid(emmap_1,emmap2,c1,c2)
      else:
        emmap1, emmap2 = match_grid(emmap1,emmap2,c1,c2)
   
    nmi = 0.0

    try:
      nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
      if nmi < 0.0: nmi = 0.0
    except:
      print('Exception for NMI score')
      print_exc()
      nmi = 0.0
    scores.append(nmi)
    idx+=1 

  return scores


def score(session, atomic_model1 = None, map_model1 = None, atomic_model2 = None, map_model2 = None, rez1 = None, rez2 = None, c1 = None, c2 = None):
  """ Generate the NMI score for 2 maps, 1 map and 1 model or 2 models. """

  sc = ScoringFunctions()
 
  # We have choices - 1 map and one model, 2 maps or 2 models
  emmap1 = None
  emmap2 = None

  blurrer = StructureBlurrer()
  
  if atomic_model1 != None and map_model1 != None:
    # 1 map 1 model
    if rez1 == None: 
      print("Please provide the resolution for the model.")
      return

    emmap1 = chimera_to_tempy_map(map_model1)
    bms = chimera_to_tempy_model(atomic_model1)
    emmap2 = blurrer.gaussian_blur(bms, rez1, densMap=emmap1)

  elif map_model1 != None and map_model2 != None:
    # 2 maps
    emmap1 = chimera_to_tempy_map(map_model1)
    emmap2 = chimera_to_tempy_map(map_model2)

  elif atomic_model1 != None and atomic_model2 != None:
    # 2 models
    if None in [rez1,rez2]: 
      print("Please provide the resolution for both model")
      return
   
    bms1 = chimera_to_tempy_model(atomic_model1)
    bms2 = chimera_to_tempy_model(atomic_model2)

    emmap1 = model_contour( bms1, rez1, emmap=False,t=0.5)
    if c2 is None: 
      emmap2 = model_contour(bms2, rez2,emmap=False,t=0.5)
    else:
      emmap2 = blur_model( bms2, rez2, emmap=False)
 
  else:
    print("Error. Must have 1 model and 1 map, 2 maps or 2 models")
    return

  # Contouring
  if c1 == None:
    c1 = map_contour(emmap1,t=1.5)

  if c2 == None:
    c2 = map_contour(emmap2,t=1.5)

  # Some kind of fix if the maps don't match?
  # Resize, resample or blur of somekind
  if not sc.mapComparison(emmap1,emmap2):
    emmap1._crop_box(c1,0.5)
    emmap2._crop_box(c2,0.5)
    
    if rez1 > 1.25*rez2: 
      emmap_2 = lpfilter(emmap2,rez1)
      emmap1, emmap2 = match_grid(emmap1,emmap_2,c1,c2)
    elif rez2 > 1.25*rez1:
      emmap_1 = lpfilter(emmap1,rez2)
      emmap1, emmap2 = match_grid(emmap_1,emmap2,c1,c2)
    else:
      emmap1, emmap2 = match_grid(emmap1,emmap2,c1,c2)
 
  nmi = 0.0

  try:
    nmi = sc.MI(emmap1,emmap2,c1,c2,1,None,None,True)
    if nmi < 0.0: nmi = 0.0
  except:
    print('Exception for NMI score')
    print_exc()
    nmi = 0.0
  return nmi
