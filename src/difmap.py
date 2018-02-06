# vim: set expandtab shiftwidth=2 softtabstop=2:

# difmap
#   "genmap" - generate a new map model
#

from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.class_arg import TempyParser

from traceback import print_exc
import os,sys
import numpy as np

from chimerax.core.map.volume import Volume
from chimerax.core.atomic.structure import AtomicStructure

from .util import chimera_to_tempy_map, tempy_to_chimera_map

from . import tool_select
from . import tool_layout


def tool_genmap(toptool):
  """ Called from the GUI. """
  # Check contours
  contour1 = None
  contour2 = None
  if toptool._widget_c1_dif.isEnabled() and toptool._widget_c1_dif.isEnabled():        
    try:
      contour1 = float(toptool._widget_c1_dif.text())
      contour2 = float(toptool._widget_c2_dif.text()) 
    except:
      toptool._show_error("Check the values fo contour1 and contour2")
      return

  # Find models
  result, map0, map1 = tool_select.select_two(toptool)

  if result:
    if isinstance(map0, Volume) and isinstance(map1, Volume):
       try:
        rez1 = float(toptool._widget_rez1_dif.text())
        rez2 = float(toptool._widget_rez2_dif.text())
        genmap(toptool.session,map0,map1,rez1,rez2,contour1,contour2)
       except:
        toptool._show_error("Check the values for rez1 and rez2")
        return
    else :
      toptool._show_error("Please select two maps.")
      return

def map_contour(m,t=-1.):
  c1 = None
  if t != -1.0:
    zeropeak,ave,sigma1 = m._peak_density()
    if not zeropeak is None: c1 = zeropeak+(t*sigma1)
    else:
      c1 = 0.0
  return c1

def genmap(session, map0 = None, map1 = None, rez1 = None, rez2 = None, c1 = None, c2 = None):
  """ Generate our new map."""
  m0 = chimera_to_tempy_map(map0)
  m1 = chimera_to_tempy_map(map1)

  # What do we do with the contours? We may already have them?
  # TODO - pull contours from m0,m1

  #MAIN CALCULATION
  #whether to shift density to positive values
  if c1 == None:
    c1 = map_contour(m0,t=1.5)
  if c2 == None:
    c2 = map_contour(m1,t=1.5)

  c1 = (c1 - m0.min())
  c2 = (c2 - m1.min())

  m0.fullMap = (m0.fullMap - m0.min())
  m1.fullMap = (m1.fullMap - m1.min())

  #find a common box to hold both maps
  spacing = max(m0.apix,m1.apix)
  grid_shape, new_ori = m0._alignment_box(m1,spacing)

  emmap_1 = m0.copy()
  emmap_2 = m1.copy()

  #resample scaled maps to the common grid
  spacing = max(rez1,rez2)*0.33

  # Not sure we should do scaling here?
  sc = ScoringFunctions()
  emmap_1.fullMap,emmap_2.fullMap = sc._amplitude_match(m0,m1,0,0,0.02,0,0,max(rez1,rez2),lpfiltb=True,lpfilta=False,ref=False)


  apix_ratio = emmap_1.apix/spacing
  diff1 = emmap_1._interpolate_to_grid(grid_shape,spacing,new_ori,1)
  diff2 = emmap_2._interpolate_to_grid(grid_shape,spacing,new_ori,1)

  # get mask inside contour for the initial maps
  emmap_1.fullMap = (m0.fullMap>c1)*1.0
  emmap_2.fullMap = (m1.fullMap>c2)*1.0
  #interpolate masks into common grid
  mask1 = emmap_1._interpolate_to_grid(grid_shape,spacing,new_ori,1,'zero')
  mask2 = emmap_2._interpolate_to_grid(grid_shape,spacing,new_ori,1,'zero')

  mask1.fullMap = mask1.fullMap > 0.1
  mask2.fullMap = mask2.fullMap > 0.1

  #min of minimums in the two scaled maps
  min1 = diff1.min()
  min2 = diff2.min()
  min_scaled_maps = min(min1,min2)
  #shift to positive values
  diff1.fullMap = diff1.fullMap - min_scaled_maps
  diff2.fullMap = diff2.fullMap - min_scaled_maps
  #range of values in the scaled maps
  min1 = np.amin(diff1.fullMap[mask1.fullMap])
  diffc1 = min1+0.10*(np.amax(diff1.fullMap)-min1)
  min2 = np.amin(diff2.fullMap[mask2.fullMap])
  diffc2 = min2+0.10*(np.amax(diff2.fullMap)-min2)
      
  #calculate difference
  diff_map = diff1.copy()
 
  #calculate difference
  diff1.fullMap = (diff1.fullMap - diff2.fullMap)
  diff2.fullMap = (diff2.fullMap - diff_map.fullMap)
 
  diff1.fullMap = diff1.fullMap*(mask1.fullMap)
  diff2.fullMap = diff2.fullMap*(mask2.fullMap)
    
  #interpolate back to original grids
  #mask1 = diff1._interpolate_to_grid1(m0.fullMap.shape,m0.apix,m0.origin,1,'zero')
  mask1 = diff1._interpolate_to_grid(m0.fullMap.shape,m0.apix,m0.origin,1,'zero')
  mask2 = diff2._interpolate_to_grid(m1.fullMap.shape,m1.apix,m1.origin,1,'zero')

  # for assigning differences (see below), use positive differences
  mask1.fullMap = mask1.fullMap*(mask1.fullMap>0.)
  mask2.fullMap = mask2.fullMap*(mask2.fullMap>0.)  

  nm0 = tempy_to_chimera_map(session, mask1)
  nm1 = tempy_to_chimera_map(session, mask2)

  session.models.add([nm0,nm1])
