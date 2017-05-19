#Added by AGNEL PRAVEEN JOSEPH
#Generate difference maps with amplitude scaling/matching, remove dusts 

import re
import sys
from TEMPy.MapParser import MapParser
from TEMPy.ScoringFunctions import ScoringFunctions
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer
#import numpy as np
import os
from TEMPy.class_arg import TempyParser
#from datetime import datetime
EXAMPLEDIR = 'Test_Files'
additional_output = False
#print help
'''
print 'use --help for help'
print '-m/-m1 [map] for input map, -m1,-m2 for two input maps'
print '-p/-p1 [pdb] for input pdb'
print '-r [resolution]; -r1,r2 for the two map resolutions'
print '-t [density threhold]; -t1,t2 for the two map thresholds'
print '--nodust to disable dusting of difference maps'
print '--softmask to apply softmask for masked map input'
print '--refscale to consider second map (or model) amplitudes as reference for scaling the first map'
print '-sw [shellwidth] to change the default shell width (1/Angstroms) for scaling (default 0.02)'
print '-dp [probability] probability of finding dust (based on size) among the difference map densities (default (>)0.2)'
print '--nofilt to disable lowpass filter before amplitude scaling (not recommended)'
print '--noscale to disable amplitude scaling (not recommended)'
'''

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
apix = tp.args.apix
#Whether to smooth unmasked region? (for hard masks)
msk = tp.args.softmask
#whether to scale amplitudes
flag_scale = True
if tp.args.noscale: 
    flag_scale = False
    print 'Warning: scaling disabled!'
#whether to lowpass filter before scaling
flag_filt = True
if tp.args.nofilt:
    flag_filt = False
#whether to use the second map (model map) as reference
refsc=tp.args.refscale
if not tp.args.mode is None: 
  if tp.args.mode == 2: refsc=True
# width of resolution shell
sw = tp.args.shellwidth
#whether to apply dust filter after difference
flag_dust = True
if tp.args.nodust: flag_dust = False
randsize = 0.2
if flag_dust:
  randsize = tp.args.dustprob
flag_softdust = False
if tp.args.softdust: flag_softdust = True
#print tp.args
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
output_synthetic_map = False
if flag_example:
    m1 = os.path.join(path_example,'emd_1046.map')
    m2 = os.path.join(path_example,'emd_1047_resampled_1046.mrc')
    r1 = 23.5
    r2 = 14.5
    Name1 = os.path.basename(m1).split('.')[0]
    Name2 = os.path.basename(m2).split('.')[0]
    c1 = 0.0607
    c2 = 0.0597
    emmap1=MapParser.readMRC(m1)
    emmap2=MapParser.readMRC(m2)
    emmap1.fullMap = (emmap1.fullMap-emmap1.mean())#emmap1.mean())/emmap1.std()
    flag_filt = False
    sw = 0.001
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
    Name2,emmap2,c2 = model_contour(p,res=r1,emmap=emmap1,t=0.5)
    #scale based on the model amplitudes
    refsc = True   
    if additional_output: output_synthetic_map = True
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
        
#MAIN CALCULATION
#whether to shift density to positive values
'''
c1 = (c1 - emmap1.min())
c2 = (c2-emmap2.min())

emmap1.fullMap = (emmap1.fullMap - emmap1.mean())
emmap2.fullMap = (emmap2.fullMap - emmap2.mean())
'''

#emmap1._crop_box(c1,2)
#emmap2._crop_box(c2,2)
#find a common box to hold both maps
spacing = max(emmap1.apix,emmap2.apix)
grid_shape, new_ori = emmap1._alignment_box(emmap2,spacing)

emmap_1 = emmap1.copy()
emmap_2 = emmap2.copy()
#if a soft mask has to be applied to both maps
if msk:
    print 'Applying soft mask'
    emmap1.fullMap = emmap1._soft_mask(c1)
    emmap2.fullMap = emmap2._soft_mask(c2)
#print datetime.now().time()

sc = ScoringFunctions()
if flag_scale:
    print 'scaling'
    if refsc: print 'Using second model/map amplitudes as reference'
    # amplitude scaling independant of the grid
    emmap_1.fullMap,emmap_2.fullMap = sc._amplitude_match(emmap1,emmap2,0,0,sw,0,0,max(r1,r2),lpfiltb=flag_filt,lpfilta=False,ref=refsc)


#resample scaled maps to the common grid
if apix is None: spacing = max(r1,r2)*0.33
else: spacing = apix
apix_ratio = emmap_1.apix/spacing
diff1 = emmap_1._interpolate_to_grid1(grid_shape,spacing,new_ori,1)
diff2 = emmap_2._interpolate_to_grid1(grid_shape,spacing,new_ori,1)

# get mask inside contour for the initial maps
emmap_1.fullMap = (emmap1.fullMap>c1)*1.0
emmap_2.fullMap = (emmap2.fullMap>c2)*1.0
mask1 = emmap_1._interpolate_to_grid1(grid_shape,spacing,new_ori,1,'zero')
mask2 = emmap_2._interpolate_to_grid1(grid_shape,spacing,new_ori,1,'zero')
mask1.fullMap = mask1.fullMap > 0.5
mask2.fullMap = mask2.fullMap > 0.5

print 'calculating difference'
# find difference map and apply contour mask
diff_map = diff1.copy()
diff1.fullMap = (diff1.fullMap - diff2.fullMap)*(mask1.fullMap)
diff2.fullMap = (diff2.fullMap - diff_map.fullMap)*(mask2.fullMap)
'''
#find level after downsample based on volume
c1 = ScoringFunctions._find_level(diff1,np.sum(emmap1.fullMap>c1)*(emmap1.apix**3))
c2 = ScoringFunctions._find_level(diff2,np.sum(emmap2.fullMap>c2)*(emmap2.apix**3))

# difference maps with contour masks
diff_map = diff1.copy()
diff1.fullMap = (diff1.fullMap - diff2.fullMap)*(diff1.fullMap>c1)
diff2.fullMap = (diff2.fullMap - diff_map.fullMap)*(diff2.fullMap>c2)
'''
#dust filter

if flag_dust:
    print 'dusting'
    #TODO: fix the contour selection for dusting
    ##diff1.fullMap = diff1._label_patches(max(diff1.min(),diff1.max()-5.0*diff1.std()),prob=randsize)[0]
    ##diff2.fullMap = diff2._label_patches(max(diff2.min(),diff2.max()-5.0*diff2.std()),prob=randsize)[0]
    if flag_softdust:
        diff1.fullMap = diff1._label_patches(diff1.min()+(diff1.max()-diff1.min())*0.5,prob=randsize)[0]
        diff2.fullMap = diff2._label_patches(diff2.min()+(diff2.max()-diff2.min())*0.5,prob=randsize)[0]
    else:
        diff1.fullMap = diff1._label_patches(diff1.min()+(diff1.max()-diff1.min())*0.75,prob=randsize)[0]
        diff2.fullMap = diff2._label_patches(diff2.min()+(diff2.max()-diff2.min())*0.75,prob=randsize)[0]
    diff1.fullMap = (diff1.fullMap>0.0)*diff1.fullMap
    diff2.fullMap = (diff2.fullMap>0.0)*diff2.fullMap

#interpolate back to original grids
mask1 = diff1._interpolate_to_grid1(emmap1.fullMap.shape,emmap1.apix,emmap1.origin,1,'zero')
mask2 = diff2._interpolate_to_grid1(emmap2.fullMap.shape,emmap2.apix,emmap2.origin,1,'zero')

'''
# get mask for difference map
diff1.fullMap = (diff1.fullMap>0.0)*1.0
diff2.fullMap = (diff2.fullMap>0.0)*1.0
mask1 = diff1._interpolate_to_grid1(emmap1.fullMap.shape,emmap1.apix,emmap1.origin,1,'zero')
mask2 = diff2._interpolate_to_grid1(emmap2.fullMap.shape,emmap2.apix,emmap2.origin,1,'zero')
#original map values if necessary # optional when the low pass filter has artifacts
mask1.fullMap = (mask1.fullMap > 0.0)*1.0
mask2.fullMap = (mask2.fullMap > 0.0)*1.0
emmap1.fullMap = emmap1.fullMap * mask1.fullMap
emmap2.fullMap = emmap2.fullMap * mask2.fullMap
'''

mask1.write_to_MRC_file(Name1+'-'+Name2+'.mrc')
mask2.write_to_MRC_file(Name2+'-'+Name1+'.mrc')

# If PDB given write out synthetic map
if output_synthetic_map:
    print 'Output synthetic map from : ', Name2
    syn_map = emmap2._interpolate_to_grid1(emmap1.fullMap.shape,
                                                emmap1.apix,
                                                emmap1.origin,
                                                1,
                                                'zero')
    syn_map.write_to_MRC_file(Name2+'_syn.mrc')

