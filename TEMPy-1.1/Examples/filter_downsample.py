import re
import sys
from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
import numpy as np
import os
from TEMPy.class_arg import TempyParser

tp = TempyParser()
tp.generate_args()

if not tp.args.inp_map is None:
  m1 = tp.args.inp_map
elif not tp.args.inp_map1 is None:
  m1 = tp.args.inp_map1
else: sys.exit('Input map (-m) and resolution (-r) required, contour (-t) optional but recommended')

if not tp.args.res is None:
  r1 = tp.args.res
elif not tp.args.res1 is None:
  r1 = tp.args.res1
else: sys.exit('Input map (-m) and resolution (-r) required, contour (-t) optional but recommended')

if not tp.args.apix is None:
  spacing = tp.args.apix
else:
  spacing = float(r1)/3.

def downsample_map(emmap,spacing):
  apix_ratio = emmap.apix/spacing
  grid_shape = int(round(emmap.x_size()*apix_ratio)),int(round(emmap.y_size()*apix_ratio)),int(round(emmap.z_size()*apix_ratio))
  try: emmap_1 = emmap._interpolate_to_grid1(grid_shape,spacing,emmap.origin)
  except: emmap_1 = emmap._interpolate_to_grid(grid_shape,spacing,emmap.origin)
  return emmap_1


print 'reading map' 
m1Name = os.path.basename(m1).split('.')[0]
emmap1=MapParser.readMRC(m1)
cutoff = emmap1.apix/float(r1)
#lowpass filter
print 'lowpass'
map_filt = emmap1._tanh_lowpass(cutoff)
#downsample after filter
print 'downsample'
emmap_1 = downsample_map(map_filt,spacing)

emmap_1.write_to_MRC_file(m1Name+'_filt_downsamp.mrc')

