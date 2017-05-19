##===============================================================================
#     This file is part of TEMPy.
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     Copyright 2010-2014 TEMPy Inventors and Birkbeck College University of London.
#                          The TEMPy Inventors are: Maya Topf, Daven Vasishtan, 
#                           Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph,
#                          Harpal Sahota
# 
# 
#     TEMPy is available under Public Licence.
#     
#     Please cite your use of TEMPy in published work:
#     
#     Vasishtan D, Topf M. (2011) J Struct Biol 174:333-343. Scoring functions for cryoEM density fitting.
#     Pandurangan AP, Vasishtan D, Topf M. (2015) Structure 23:2365-2376. GAMMA-TEMPy: Simultaneous fitting of components in 3D-EM Maps of their assembly using genetic algorithm.
#===============================================================================

#!/Users/arun/Library/Enthought/Canopy_64bit/User/bin/python
from TEMPy.StructureParser import PDBParser
from TEMPy.MapParser import MapParser
#from TEMPy.StructureBlurrer import StructureBlurrer
#from TEMPy.ScoringFunctions import *
#from TEMPy.Assembly import Assembly
from TEMPy.Gamma import *
from TEMPy.VQ import *
import time
import getopt

def usage():
	print "\nGAMMA-TEMPy version 1.1\n"
	print "Usage:\n"
	print 'multiple_ga --ipdb <pdbfile> --imap <mrcfile> --res <resolution> [--options]'
	print "--popsize Number of members in the population (default=160)\n"
	print "--ngen Number of generations in the GA (default=100)\n"
	print "--nga Number of GA solutions to generate (default=1)\n"
	print "--gof Goodness-of-fit score to use. (1)Mutual information score, (2)CCC (default=1)\n"
	print "--ncpu Number of CPU to use (default=1)\n"
	print "--outdir Name of the output directory to store results (default is set to the name of the input PDB file name)\n"
	print "--outfile Prefix name of the predicted fits (default is set to the name of the input PDB file name)\n"
	print "--pdbsnap (1) To write all the members (PDB fits) in every GA population for every generation. (default=0; Do not write PDB fits in the population)\n"
	print "--seedpoints Supply a PDB formated file containing the coordinate to use as initial VQ points (By default the VQ points will be generated from the input density map)\n"
	print "--mrate Mutation rate (default=0.02). This is a advanced setting for expert user\n"
	print "--crate Crossover rate (default=0.8). This is a advanced setting for expert user\n"
	print "--moffset Magnitude of the displacement used while applying transational mutation (default set 0.0 will calculate moffset value from the VQ point distances). This is a advanced setting for expert user\n"
	print "NOTE: Parallel processin is done using Parallel Python (PP). If more machines and nodes are required, then start the ppserver in auto discovery mode (-a -d option) in the remote nodes\n"
	print "See www.parallepython.com for further information\n"


input_config = {"--help":'', "--ipdb":'dummy.pdb', "--imap":'dummy.mrc', "--res":0, "--popsize":160, "--ngen":100, "--nga":1, "--gof":1, "--ncpu":1, "--outdir":'dummy', "--outfile":'dummy.log', "--pdbsnap":0, "--seedpoints":'null', "--mrate":0.02, "--crate":0.8, "--moffset":0.0} 
try:
	opts, args = getopt.getopt(sys.argv[1:],'',["ipdb=","imap=","res=","popsize=","ngen=","nga=","gof=","ncpu=","outdir=","outfile=","pdbsnap=","seedpoints=","mrate=", "crate=", "moffset="])
	#opts, remainder = getopt.getopt(sys.argv[1:],'i:o:',['ipdb=','omap='])
	input_config.update(dict(opts))
	#print input_config
	#print opts, args

except getopt.GetoptError, err:
	usage()
	#print 'multiple_ga --ipdb <pdbfile> --imap <mrcfile> --res <resolution> --popsize <popsize> --ngen <number_of_generations> --nga <number_of_ga_runs> --gof <score_type> --ncpu <number_of_cpus> --outdir <directory_name> --outfile <prefix_filename> --pdbsnap <1 or 0>' 

if input_config["--help"] == '':
	usage()

if input_config["--ipdb"] == "dummy.pdb":
	print 'Input pdb name not specified.\n'
	usage()
	#print 'multiple_ga --ipdb <pdbfile> --imap <mrcfile> --res <resolution> --popsize <popsize> --ngen <number_of_generations> --nga <number_of_ga_runs> --gof <score_type> --ncpu <number_of_cpus> --outdir <directory_name> --outfile <prefix_filename> --pdbsnap <1 or 0>' 
	exit(0)

if input_config["--imap"] == "dummy.mrc":
	print 'Input map name not specified.\n'
	usage()
	#print 'multiple_ga --ipdb <pdbfile> --imap <mrcfile> --res <resolution> --popsize <popsize> --ngen <number_of_generations> --nga <number_of_ga_runs> --gof <score_type> --ncpu <number_of_cpus> --outdir <directory_name> --outfile <prefix_filename> --pdbsnap <1 or 0>' 
	exit(0)

if input_config["--res"] == 0:
	print 'Input resolution value not specified.\n'
	usage()
	#print 'multiple_ga --ipdb <pdbfile> --imap <mrcfile> --res <resolution> --popsize <popsize> --ngen <number_of_generations> --nga <number_of_ga_runs> --gof <score_type> --ncpu <number_of_cpus> --outdir <directory_name> --outfile <prefix_filename> --pdbsnap <1 or 0>' 
	exit(0)

if input_config["--outdir"] == 'dummy':
	print input_config["--ipdb"]
	input_config["--outdir"] = input_config["--ipdb"].split('.')[0]
if input_config["--outfile"] == 'dummy.log':
	input_config["--outfile"] = input_config["--ipdb"].split('.')[0]+'.log'
	
# Read pdb
try:
	prot = PDBParser.read_PDB_file(input_config["--outdir"],input_config["--ipdb"],hetatm=False,water=False)
	# Get number of components
	comps = prot.split_into_chains()
	ncomp = len(comps)
except:
	print "Error in reading pdb coordinate file\n"
	usage()
	exit(0)
	
try:
	emmap = MapParser.readMRC(input_config["--imap"])
	emmap.normalise()
except:
	print "Error in reading map density file\n"
	usage()
	exit(0)
try:
	res = float(input_config["--res"])
	if res == 0.0:
		raise
except:
	print "Error in reading map resolution value\n"
	usage()
	exit(0)


#Set the machine names for prallel processing
#Please do not change the below tuple for the purpose of the practicle
#machines = ("talos","emma","ceres","eureka")
#machines = ("talos",)

# Initialise GA
ga = GA()
total_popsize = int(input_config["--popsize"]) # Set population size
ngen = int(input_config["--ngen"]) # Set number of GA generations
grun = int(input_config["--nga"]) # Set number of independent GA runs
selection_method = 1 # Set selection method; 1=Tournament selection
gof = int(input_config["--gof"]) # Choose the goodness-of-fit score to use; 1= Mutual information, 2= CCC
w_gof = ncomp #Set the weight for gof score
w_clash = 1.0 #Set the weight for clash score
mrate = float(input_config["--mrate"])
crate = float(input_config["--crate"])
moffset = float(input_config["--moffset"])
ncpu = int(input_config["--ncpu"]) # Set the number of cores to use in parallel

# Set the prefix for the output filename
output_dir = input_config["--outdir"]
output_filename = input_config["--outfile"]

if os.path.isdir(output_dir):
	print 'The output directory already exist. Please give a different directory name'
	exit(0)
else:
	os.mkdir(output_dir)

output_name = output_dir+"/"+output_filename	

# Set the pdb filename prefix to write out the fit from the GA generations

#pdb_snapshot_dir = output_dir+'_snapshot'
#pdb_snapshot_file = input_config["--outfile"]

#try:
#	os.stat(pdb_snapshot_dir)
#except:
#	os.mkdir(pdb_snapshot_dir)
#pdb_snapshot_name = pdb_snapshot_dir+"/"+pdb_snapshot_file

try:
	if int(input_config["--pdbsnap"]) != 0:
		pdb_snapshot_dir = output_dir+'_snapshot'
		pdb_snapshot_file = input_config["--outfile"]

		try:
			os.stat(pdb_snapshot_dir)
		except:
			os.mkdir(pdb_snapshot_dir)
		pdb_snapshot_name = pdb_snapshot_dir+"/"+pdb_snapshot_file
	else:
		pdb_snapshot_name = 'dummy'
except:
	print 'Error in pdbsnap parameter.\n'
	usage()
	exit(0)

if input_config["--seedpoints"] == 'null':
	#print 'taking seed points from map'
	vq_vec_list = get_VQ_points(emmap, emmap.std()*2, ncomp, 50, output_name, False)
else:
	vq_struct = PDBParser.read_PDB_file('vq',input_config["--seedpoints"],hetatm=False,water=False)
	vq_vec_list = vq_struct.get_vector_list()
	#print 'taking seedpoints from pdb'
	if len(vq_vec_list) < ncomp :
		print 'Not enough seed coordinates supplied. Ending program.\n'
		print "Printing the VQ point list\n"
		print vq_vec_list
		exit(0)

f = file(output_name+'.log', 'a')
start_time = time.time()
f.write("Start time : "+str(start_time)+"\n")
f.write("------------------------------------------------------------------------------------------\n")
f.write("GAMMA-TEMPy\n")
f.write("Genetic Algorithm for Modelling Macromolecular Assemblies\n")
f.write("------------------------------------------------------------------------------------------\n")
f.write("Input assembly                : ../"+str(input_config["--ipdb"])+"\n")
f.write("Input map                     : ../"+str(input_config["--imap"])+"\n")
f.write("Resolution (angstrom)         : "+str(input_config["--res"])+"\n")
#f.write("Population size               : "+str(total_popsize)+"\n")
#f.write("Number of GA Generations      : "+str(ngen)+"\n")
#f.write("Number of GA runs             : "+str(grun)+"\n")
f.write("Mutation rate                 : "+str(mrate)+"\n")
f.write("Crossover rate                : "+str(crate)+"\n")
f.write("Mutation ofset                : "+str(moffset)+"\n")
if gof == 1:
	f.write("Goodness-of-fit score         : Mutual Information\n")
else:
	f.write("Goodness-of-fit score         : Cross Correlation Coefficient\n")
#f.write("Output directory : "+output_dir+"\n")
f.write("Prefix of the output filename : "+output_filename+"\n")
if int(input_config["--pdbsnap"]) != 0:
	f.write("Output assembly model after every generation : Yes ("+pdb_snapshot_name+")\n")
else:
	f.write("Output assembly model after every generation : No\n")

if input_config["--seedpoints"] == 'null':
	f.write("Map vector quantisation file  : "+output_filename+'_vq.pdb'+"\n")
else:
	f.write("Map vector quantisation file  : "+input_config["--seedpoints"]+"\n")
#f.write("Number of CPUs used for the calculation : "+str(ncpu)+"\n")
f.close()
#Run GA assembly fitting
#results = ga.run(grun, ngen, total_popsize, selection_method, gof, w_gof, w_clash, prot, ncomp, emmap, res, output_name, pdb_snapshot_name, vq_vec_list, machines, mrate, crate, moffset, ncpu)
results = ga.run(grun, ngen, total_popsize, selection_method, gof, w_gof, w_clash, prot, ncomp, emmap, res, output_name, pdb_snapshot_name, vq_vec_list, mrate, crate, moffset, ncpu)

#print 'Seconds : ', time.time()-start_time
f = file(output_name+'.log', 'a')
f.write("------------------------------------------------------------------------------------------\n")
f.write("Execution time (sec): "+str(time.time()-start_time)+"\n")
f.close()
exit(0)


