import re
import sys
from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
#import numpy as np
import os
from TEMPy.class_arg import TempyParser

tp = TempyParser()
tp.generate_args()

if not tp.args.inp_map is None:
  m1 = tp.args.inp_map
elif not tp.args.inp_map1 is None:
  m1 = tp.args.inp_map1
else: sys.exit('Input map missing')
print 'reading map' 
m1Name = os.path.basename(m1).split('.')[0]
emmap_1=MapParser.readMRC(m1)

if not tp.args.thr is None: c1 = tp.args.thr
elif not tp.args.thr1 is None: c1 = tp.args.thr1
else:
  print 'calculating contour'
  zeropeak,ave,sigma1 = emmap_1._peak_density()
  if not zeropeak is None: c1 = zeropeak+(1.5*sigma1)
  else:
    sys.exit('Contour level required')


level = c1
sigma = emmap_1.fullMap.std()
sigma = abs(sigma)
try: emmap_1.fullMap = emmap_1._label_patches(level-0.02*sigma)[0]
except: emmap_1._map_binary_opening(level-0.02*sigma)
##emmap_1._mask_contour(level,0.5)
emmap_1.write_to_MRC_file(m1Name+'_contour_dust.mrc')

