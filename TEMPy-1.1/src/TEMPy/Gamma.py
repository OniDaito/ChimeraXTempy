##===============================================================================
#     This file is part of TEMPy.
#     It implements a genetic algorithm for the purpose of fitting multiple component into the assembly map
#     
#     TEMPy is a software designed to help the user in the manipulation 
#     and analyses of macromolecular assemblies using 3D electron microscopy maps. 
#     
#     This module handles the fitting of multiple subunits into EM maps using Genetic Algorithm	
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

from random import randint, randrange, uniform, sample
from TEMPy.Vector import *
from copy import deepcopy
import cPickle as pickle
from TEMPy.Quaternion import *
from TEMPy.ScoringFunctions import *
from TEMPy.StructureBlurrer import StructureBlurrer
from TEMPy.StructureParser import PDBParser
from TEMPy.Assembly import *
from TEMPy.VQ import *
import numpy as np
import pp
import cPickle as pickle
import os
from time import *
import subprocess
import multiprocessing

# Get the euler angle dictionary from the pickle file
this_dir, this_filename = os.path.split(__file__)
#print this_dir, this_filename
DATA_PATH = os.path.join(this_dir, "tempy_data", "quaternion_vectors_5000.pk")

quat_vec = pickle.load(open(DATA_PATH, "r"))

def move_struct_quat(g,assembly):
    """ 
    
    Method used to apply the translation and rotation opetations on the individual componets using the information in a given genotype.

	Arguments:
	    *g*
		An instance of a Genotype object.
	    *assembly*
		An instance of the Assembly object. 
    
    """

    from copy import deepcopy
    for vec in range(0,len(g),2):
	tx = g[vec].x 
	ty = g[vec].y
	tz = g[vec].z
	q_param = g[vec+1].param
	mat = g[vec+1].to_rotation_matrix()
	assembly.move_map_and_prot_by_quat(vec/2, tx, ty, tz, q_param, mat)
    newstruct = deepcopy(assembly.combine_structs())
    newmap = deepcopy(assembly.combine_maps())
    newmaplist = deepcopy(assembly.mapList)
    assembly.reset_all()
    return newmap, newstruct, newmaplist

def score_pop_segment_MI(pop_seg, scorer, assembly, emmap, ncomp, cvol, template_grid, apix, w_gof, w_clash):
    """ 
    
    Method used to score the genotypes (fits) in the population.

	Arguments:
	    *pop_seg*
	        An instance of a Population object.
	    *scorer*
	        An instance of the ScoringFunction object. 
	    *assembly*
	        An instance of the Assembly object. 
	    *emmap*
	        An instance of the map object.
	    *ncomp*
	        Number of components in the assembly.
 	    *cvol*
		List containing the volume values of the individual components. Used in the calculation of the clash score.
	    *template_grid*
		Map instance to be used as a template for the purpose of calculating the clash score.
	    *apix*
		voxel size of the grid
	    *w_gof*
	        Weighting used for Goodness-of-fit score contribution to the GA fitness score. 
	    *w_clash*
	        Weighting used for clash penalty score contribution to the GA fitness score.

	Return:
	    Return an instance of a Population object with all its member's fitness been score.         
   
    """

    from itertools import combinations
    for g in pop_seg:
	gof = 0.0
	clash_score = 0.0
	protruction = 0
	if g.fitness == 0:
		newmap, newstruct, newmaplist = move_struct_quat(g,assembly)		
		gof = scorer.MI(newmap, emmap)*w_gof
        	clash_score = scorer.get_sm_score(newstruct, ncomp, template_grid, cvol, apix)*w_clash
        	symscore = 0.0
        	g.gof += gof
        	g.clash += clash_score
        	g.protruction += protruction
        	g.symres += symscore
        	g.fitness += gof + clash_score + protruction + symscore
    return pop_seg

def score_pop_segment_CCC(pop_seg, scorer, assembly, emmap, ncomp, cvol, template_grid, apix):
    """ 
    
    Method used to score the genotypes (fits) in the population.

	Arguments:
	    *pop_seg*
	        An instance of a Population object.
	    *scorer*
	        An instance of the ScoringFunction object. 
	    *assembly*
	        An instance of the Assembly object. 
	    *emmap*
	        An instance of the map object.
	    *ncomp*
	        Number of components in the assembly.
 	    *cvol*
		List containing the volume values of the individual components. Used in the calculation of the clash score.
	    *template_grid*
		Map instance to be used as a template for the purpose of calculating the clash score.
	    *apix*
		voxel size of the grid
	    *w_gof*
	        Weighting used for Goodness-of-fit score contribution to the GA fitness score. 
	    *w_clash*
	        Weighting used for clash penalty score contribution to the GA fitness score.

	Return:
	    Return an instance of a Population object with all its member's fitness been score.         
   
    """

    from itertools import combinations
    for g in pop_seg:
	gof = 0.0
	clash_score = 0.0
	protruction = 0
	if g.fitness == 0:
		newmap, newstruct, newmaplist = move_struct_quat(g,assembly)		
		gof = scorer.CCC(newmap, emmap)*w_gof
        	clash_score = scorer.get_sm_score(newstruct, ncomp, template_grid, cvol, apix)*w_clash
        	symscore = 0.0
        	g.gof += gof
        	g.clash += clash_score
        	g.protruction += protruction
        	g.symres += symscore
        	g.fitness += gof + clash_score + protruction + symscore
    return pop_seg

class GA():
        """ 
    
        A class used for implementing the Genetic Algorithm (GA). 
    
        """

	def __init__(self):
		pass

	#def run(self, runs, no_of_gen, pop_size, selection_method, gof, w_gof, w_clash, prot, ncomp, emmap, resolution, logfile, gasnap, vq_vec_list, machines, mrate, crate, moffset, ncpu=0):
	def run(self, runs, no_of_gen, pop_size, selection_method, gof, w_gof, w_clash, prot, ncomp, emmap, resolution, logfile, gasnap, vq_vec_list, mrate, crate, moffset, ncpu=1):
		"""
		
		Main method to initiate GA cycle.

		Arguments:
		    *runs*
		        Number of GA solution to generate.
		    *no_of_gen*
		        Number of GA generations to generate. 
		    *pop_size*
		        Number of members in the GA population. 
		    *selection_method*
		        Selection method used to pick members in the population for the purpose of generating new population.
		        Currently should be set to 1  for tournament selection.
 		    *gof*
		        Option to specify the Goodness-of-fit function to use.
			Set it to 1 for Mutual information score or 2 for Cross Correlation Coefficient score. 
		    *w_gof*
		        Weighting used for Goodness-of-fit score contribution to the GA fitness score. 
		    *w_clash*
		        Weighting used for clash penalty score contribution to the GA fitness score. 
		    *prot*
		        Instance of a Structure_BioPy object that contain multiple chains used as an input for building Assembly object 
		    *ncomp*
		        Number of component in the assembly. 
		    *emmap*
		        Instance of a map object.
		    *resolution*
		        Resolution of the map.
		    *logfile*
		        Name of the output logfile.
		    *gasnap*
		        Option used to control the PDB files written during the GA run.
			Set it to 1 for writing each individual member in the population (fit) in the every GA generation.
			Default is set to 'dummy' which will not write each individual member in the population in every GA generation.
 		    *vq_vec_list*
		        List of Vector objects used to represent initial point configuration which is used to generate initial populatons of fits.
 		    *mrate*
		        Mutation rate for the mutation operator.
 		    *crate*
		        Crossover rate for the mutation operator.
 		    *moffset*
		        Set the translation offset range (in Angstrom) applied to each of the components position generated in the initial population pool.
 		    *ncpu*
		        Number of cpus to use in parallel through Parallel Python.
		Return:
		    The function return the following items.
		    An instance of the Population object of the final generation.
	            A Structure_BioPy object corresponding to the fittest member in the final genaration and its respective simulated map object.
	            A string containing the best fitness score, Min, Max, Avg, Std and total fitness score of all fits in the final genetation.         
		    
		"""


		sel = Selection()
		scorer = ScoringFunctions()

		if selection_method == 1:
			sel_type = sel.tournament
			tour_size = 2
		else:
			print 'Selection method other than tournament selection is not tested. Please use tournament selection'
			exit(0)

		if gof != 1 and gof != 2:
			print 'Please select the gof=1 for MI or gof=2 for CCC'
			exit(0)
			
		#Build the assembly
		assembly = Assembly(prot.split_into_chains())

		# Building maps based on emmap
		assembly.build_maps(resolution,emmap)

		
		#Template grid for volume clash score
		apix = 3.5 # apix used to record the atom position overlayed in the template map
		
		template_grid = emmap._make_clash_map(apix)

		#for grid based clash score
		#template_grid = emmap.make_clash_cplx_map(apix)
		
		#Create a blurrer object
		blurrer = StructureBlurrer()

		# Get the volumes occupied by the components in the assembly
		cvol = emmap._get_component_volumes(prot, apix, blurrer)

		#ppservers=("emma","talos")
		#If required more machines, start the ppserver in auto discovery mode (-a -d option) in the remote nodes
		ppservers = ("*",)
		
		'''
		if ncpu == 0 or ncpu > 5:
			#For running in the local machine
			#jobserver = pp.Server(ppservers=ppservers)
			#cpu_avil = multiprocessing.cpu_count()
			#For the purpose of EMBO practical
			print 'ncpu should not be 0 and should not be > 5. Exiting'
			exit(0)			
		else:
			#For running with using n cpus (a general case where we have got more remote processing machines)
			jobserver = pp.Server(ncpus=ncpu, ppservers=ppservers)
			cpu_avil = ncpu
		'''
		jobserver = pp.Server(ncpus=ncpu, ppservers=ppservers)
		cpu_avil = ncpu

		if pop_size <= cpu_avil:
			cpu_avil = pop_size
		else:
			#Used to round off the population size 
			#and make it divisible by the available number of processor
			#if ncpu != 1:
			#	pop_size = (pop_size/cpu_avil)*cpu_avil + cpu_avil
			if pop_size%cpu_avil != 0:
				pop_size = (pop_size/cpu_avil)*cpu_avil + cpu_avil
				
		if pop_size < 10:
			print 'Populations size < 10. Please increase the size.'
			exit(0)
		
		#Send the dynamically determined input parameters to log file
		f = file(logfile+'.log', 'a')
		f.write("Population size               : "+str(pop_size)+"\n")
		f.write("Number of GA Generations      : "+str(no_of_gen)+"\n")
		f.write("Number of GA runs             : "+str(runs)+"\n")
		f.write("Number of CPUs used for the calculation : "+str(cpu_avil)+"\n")
		f.write("------------------------------------------------------------------------------------------\n")
		f.write("GA_run and Generation no., Best fitness score, Weighted MI, Clash penalty, Protrusion penalty, Symmetry score, Worse fitness score, Total population fitness score, Average population fitness score, Std of population fitness score \n")
		f.close()
		
		#Testing for automatice remote logging and starting a ppserver
		'''
		secret = 'HAI'
		hosts = list(machines) #["talos"]
		ppservers = machines #("talos",)

		serv_arg1 = '/d/d610/s/epd/pp-1.6.1/ppserver.py -w 10 -s ' + secret + ' -t 300 < /dev/null > /dev/null &'
		#serv_arg1 = '/d/d610/s/epd/pp-1.6.1/ppserver.py -t 300 &'
		

		for s in hosts:
			os.spawnv(os.P_WAIT, '/usr/bin/ssh',['',s,'',serv_arg1])
			#ssh  = subprocess.Popen(["/usr/bin/ssh", "%s" % s, serv_arg1])

		print 'ppservers started'
		#sleep(15)
		jobserver = pp.Server(ncpus=0, ppservers=ppservers, secret=secret)
		#jobserver = pp.Server(ncpus=0, ppservers=ppservers)

		print jobserver.get_ncpus()
		print jobserver.get_active_nodes()
		#exit(0)
		'''
		
		# Build maps of subunits, if not already done
		if not assembly.mapList:
		    assembly.build_maps(resolution, emmap)

		n_vq = len(vq_vec_list)
		
		if moffset == 0.0:
			#Find the minimum distance among all the vq pairs
			#Set the minimum distance to max_trans that is applied around each vq point to generate random fit
			dist_list = np.zeros(sum(range(1,len(vq_vec_list))))
			p = 0
			for i in range(len(vq_vec_list)-1):
		        	for j in range(i+1, len(vq_vec_list)):
		        	        v1 = vq_vec_list[i]
		        	        v2 = vq_vec_list[j]
		        	        dist_list[p] = v1.dist(v2)
		        	        p = p + 1
			max_change = dist_list.min()
		else:
			max_change = moffset
						
		#GA run loop
		for grun in xrange(runs):
			# Get the ga_pool
			pop = self.get_ga_pool(assembly, pop_size, max_change, vq_vec_list, quat_vec)
			curr_pop = self.score_population(pop, pop_size, gof, w_gof, w_clash, cpu_avil, jobserver, scorer, assembly, emmap, ncomp, cvol, template_grid, apix)

			# Start iterations
			for x in xrange(no_of_gen):
				f = file(logfile+'.log', 'a')
				#Mutation rate linearly decreasing from 0.2 to 0.02
				mutRate = mrate+(1-x/float(no_of_gen))*0.08
				crossRate = crate

				#breed new population
				new_pop = curr_pop.breedNewPop(no_of_gen, x, mutRate, crossRate, sel_type, tour_size)
				pop_size = new_pop.size()
				new_pop = self.score_population(new_pop, pop_size, gof, w_gof, w_clash, cpu_avil, jobserver, scorer, assembly, emmap, ncomp, cvol, template_grid, apix)
				old_new_Pop = Population()
				for pold in curr_pop.pop:
					old_new_Pop.addGenotype(pold.copy())
				for pnew in new_pop.pop:
					old_new_Pop.addGenotype(pnew.copy())
	
				best = old_new_Pop.pickSetOfBest(pop_size)

				curr_pop = Population()        
				for b in best:
					curr_pop.addGenotype(b.copy())

				#Add info to log file to record the generation number,
				#best pop fittness values (total, MI, clash), min_pop_fitness, total_pop_fitness
				#average_pop_fitness and std_pop_fitness
				f.write("R"+str(grun+1)+"G"+str(x+1)+","+str(curr_pop.getBestScores())+","+str(curr_pop.min_fitness())+","+str(curr_pop.totalFitnessScore())+'\n')
				f.write(str(curr_pop)+'\n\n')

			 	if gasnap != 'dummy':
					#Writing out assembly models for every generation
					tmap, tstruct, maplist = move_struct_quat(curr_pop.pickBest(),assembly)
					tstruct.write_to_PDB(gasnap+'_'+str(grun+1)+'_'+str(x+1)+'.pdb')

			newmap, newstruct, maplist = move_struct_quat(curr_pop.pickBest(),assembly)		 	
			newstruct.write_to_PDB(logfile+'_'+str(grun+1)+'.pdb')
		f.close()
		jobserver.destroy()
		#Return final population, best fit, map of best fit
		return curr_pop, newstruct, newmap, "Generation"+str(x)+": "+str(curr_pop.getBestScores())+", "+str(curr_pop.min_fitness())+", "+str(curr_pop.totalFitnessScore())+'\n'

	def score_population(self, pop, pop_size, gof, w_gof, w_clash, cpu_avil, jobserver, scorer, assembly, emmap, ncomp, cvol, template_grid, apix):
		"""
		
		Method used to score memeber in the population.

		Arguments:
		    *pop*
		        Instance of the Population object.
		    *pop_size*
		        Number of members in the population. 
		    *gof*
		        Option to specify the Goodness-of-fit function to use.
			Set it to 1 for Mutual information score or 2 for Cross Correlation Coefficient score. 
		    *w_gof*
		        Selection method used to pick members in the population for the purpose of generating new population.
		        Currently should be set to 1  for tournament selection.
 		    *gof*
		        Option to specify the Goodness-of-fit function to use.
			Set it to 1 for Mutual information score or 2 for Cross Correlation Coefficient score. 
		    *w_gof*
		        Weighting used for Goodness-of-fit score contribution to the GA fitness score. 
		    *w_clash*
		        Weighting used for clash penalty score contribution to the GA fitness score. 
		    *cpu_avil*
		        Number of cpus to use in parallel through Parallel Python.
		    *jobserver*
		        Instance of the Jobserver object used by Parallel Python.
		    *scorer*
		        Instance of the ScoringFunction object.
		    *assembly*
		        Instance of a Assembly object.
		    *emmap*
		        Instance of a map object.
		    *ncomp*
		        Number of component in the assembly.
		    *cvol*
 			List containing the volume values of the individual components. Used in the calculation of the clash score.
		    *template_grid*
 			Map instance to be used as a template for the purpose of calculating the clash score.
		    *apix*
			voxel size of the grid
		Return:
		    Return an instance of a Population object with all its member's fitness been score.         


		"""

		#Score the populations using parallel processing
		seg_size = pop_size/cpu_avil
		
		segs = []
		for n in range(cpu_avil-1):
			segs.append(pop[n*seg_size:(n+1)*seg_size])
		segs.append(pop[(cpu_avil-1)*seg_size:])

		jobs = []
		funclist = (move_struct_quat,)
		#modlist = ("TEMPy","numpy","copy","pp","cPickle","time","random","os")
		modlist = ()
		
		if gof == 1:
			for s in segs:
				jobs.append(jobserver.submit(score_pop_segment_MI, (s, scorer, assembly, emmap, ncomp, cvol, template_grid, apix, w_gof, w_clash), funclist))
		else:
			for s in segs:
				jobs.append(jobserver.submit(score_pop_segment_CCC, (s, scorer, assembly, emmap, ncomp, cvol, template_grid, apix, w_gof, w_clash), funclist))
			
		jobserver.wait()
		scored_pop = Population()

		for j in jobs:
			scored_pop.pop.extend(j())
		return scored_pop.copy()



	def get_ga_pool(self, assembly, pop_size, max_trans, vq_vec_list, quat_vec):
		"""
		
		Method used to generate initial population for running GA.

		Arguments:
		    *assembly*
		        Instance of a Assembly object.
		    *pop_size*
		        Number of members in the population.
		    *max_trans*
		        Set the translation offset range (in Angstrom) applied to each of the components position generated in the initial population pool.
		    *vq_vec_list*
		        List of Vector objects used to represent initial point configuration which is used to generate initial populatons of fits.
		    *quat_vec*
		        List of Quaternion vectors (in the form of [w,x,y,z]) used for the purpose of applying rotation operation to the components in the assembly.

		Return:
		    Return an instance of a Population object.         


		"""

		# 50% of the members are produced by randomly placing the subunits around the vq points
		# 50% of the members are produced by randomly placing the subunits around the vq points and translating within (+max_tans1,-max_trans1)
		pop_size1 = int(pop_size*0.50)
		pop_size2 = int(pop_size*0.50)
		n_vq = len(vq_vec_list)
		nsubunit = len(assembly.structList)
		complex_list = []
		scores = []
		dx = 0
		dy = 0
		dz = 0

		for i in range(0,pop_size1):
			subcomplex_list = []
			rlist = sample(range(n_vq),nsubunit)
			rand_subunit_indx = randint(0,nsubunit)
			for j in range(0,nsubunit):
				offset = vq_vec_list[rlist[j]] - assembly.structList[j].CoM
		        	tx = offset.x
		        	ty = offset.y
		        	tz = offset.z
		        	r_index = randint(1, 5000)
		        	q_param = quat_vec[r_index]
		     		trans_list = [tx,ty,tz]
		     		quat_list = q_param
		     		subcomplex_list.append([trans_list,quat_list])
			complex_list.append(subcomplex_list)

		trans = max_trans
		for i in range(0,pop_size2):
			subcomplex_list = []
			rlist = sample(range(n_vq),nsubunit)
			rand_subunit_indx = randint(0,nsubunit)
			for j in range(0,nsubunit):
				offset = vq_vec_list[rlist[j]] - assembly.structList[j].CoM
				dx = uniform(-trans,trans)
		    		dy = uniform(-trans,trans)
		        	dz = uniform(-trans,trans)
		        	tx = offset.x + dx
		        	ty = offset.y + dy
		        	tz = offset.z + dz
		        	r_index = randint(1, 5000)
		        	q_param = quat_vec[r_index]
		        	trans_list = [tx,ty,tz]
		        	quat_list = q_param
		        	subcomplex_list.append([trans_list,quat_list])
			complex_list.append(subcomplex_list)

		# Build new population
		pop = Population()

		#construct cartesian and quaternion gene
		g_type_trans = []
		g_type_trans.append('cartesian')
		g_type_trans.append([-max_trans,max_trans])
		g_type_ri = []
		g_type_ri.append('quaternion')
		g_type_ri.append([1,5000])
		
		for x in sample(range(len(complex_list)),pop_size):
			assembly_gene = complex_list[x]
			geno = []
			trans = []
			quat = []
			for j in assembly_gene:
				trans.extend(g_type_trans)
				quat.extend(g_type_ri)
				trans.append(j[0])
				quat.append(j[1])
				geno.append(trans)
				geno.append(quat)
				trans = []
				quat = []
			pop.addGenotype(Genotype(geno))

		return pop

class VectorGene:
    """ 
    
    A class used to record the three dimensional position of a component in the assembly and to apply mutation and crossover genetic operators on them. 
    
    """


    def __init__(self, ul_list, x, y, z):
	"""
	
	Constructor to initialise the Vector gene.

	Arguments:
	    *ul_list*
	        A list containing the min and max translation allowed for the VectorGene.
	    *x*
	        x coordinate of a component.
	    *y*
	        y coordinate of a component.
	    *z*
	        z coordinate of a component.

	"""

	self.ul_list = ul_list
	self.x = x
	self.y = y
	self.z = z
	self.value = Vector(x,y,z)

    def get_gene_list(self):
	"""

	Returns a list of x,y,z coordinates of the vector gene.

	"""

	glist = list()
	glist.append(self.x)
	glist.append(self.y)
	glist.append(self.z)
	return glist

    def copy(self):
	"""

	Returns a copy of a VectorGene object.

	"""

        gene = VectorGene(self.ul_list,self.x,self.y,self.z)
        return gene
    	
    def mutate(self):
	"""

	Mutation operator that modify the x, y, and z coordinates of the VectorGene.

	"""

	r = random_vector2(self.ul_list)
       	self.value += r
       	self.x = self.value.x
       	self.y = self.value.y
       	self.z = self.value.z
       	

    def check_for_mutations(self, mutationRate):
	"""

	Method to apply mutation operator on VectorGene based on mutationRate.

	Arguments:
	    *mutationRate*
	        Rate of mutation used to decide an application of the mutation operation.

	"""

        if uniform(0,1) < mutationRate:
            self.mutate()            

    def crossover(self, otherGene, crossRate):
	"""

	Method to apply crossover operator on VectorGene based on crossover rate.

	Arguments:
	    *crossRate*
	        Rate of crossover used to decide an application of the crossover operation.

	"""

	r = uniform(0,1)
        if r <= crossRate:
	    newgene = self.copy()
	    newgene.value = otherGene.value
	    newgene.x = newgene.value.x
	    newgene.y = newgene.value.y
	    newgene.z = newgene.value.z
            return newgene.copy()
        else:
            return self.copy()

    def __repr__(self):
	"""

	Returns the 3D coordinate of the VectorGene as a Vector object.

	"""

        return str(self.value)

class QuaternionGene:
    """ 
    
    A class used to represent the rotational state a component in the assembly and to apply mutation and crossover genetic operators on the component. 
    
    """

    def __init__(self, ul_list, w, x, y, z):
	"""
	
	Constructor to initialise the Quaternion gene.

	Arguments:
	    *ul_list*
	        A list containing the start and the end index of the uniformly distributed quaternion.
	    *w*
	        Rotation angle component of the QuaternionGene.
	    *x*
	        x axis of the quaternion.
	    *y*
	        y axis of the quaternion.
	    *z*
	        z axis of the quaternion.

	"""

	self.ul_list = ul_list
	self.param = [w,x,y,z]

    def __repr__(self):
	"""

	Returns the components of the quaternion objects as a list object.

	"""

	return str(self.param)
 
    def copy(self):
	"""

	Returns a copy of a QuaternionGene object.

	"""

	gene = QuaternionGene(self.ul_list, self.param[0], self.param[1], self.param[2], self.param[3])
	return gene

    def crossover(self, otherquat, crossRate):
	"""

	Method returns a QuaternionGene object among two QuaternionGene based on crossover rate.

	Arguments:
	    *otherquat*
	        An instance of a QuaternionGene object.

	"""

	r1 = uniform(0,1)
	if r1 <= crossRate:
	    return otherquat.copy()
        else:
	    return self.copy()
	
    def get_interpolated_quat(self, q2_param):
	"""

	Return an interpolated QuaternionGene found between two different QuaternionGene.

	Arguments:
	    *q2_param*
	        A list of type [w,x,y,z] used of the purpose of finding an interpolated quaternion between two quaternions.

	"""

	w1, x1, y1, z1 = self.param
        w2, x2, y2, z2 = q2_param

        fraction = uniform(0,1)
	lamda = self.dot_product(q2_param)
        if lamda < 0:
                w2 = -w2
                x2 = -x2
                y2 = -y2
                z2 = -z2
        w = w1 + fraction * (w2 - w1)
        x = x1 + fraction * (x2 - x1)
        y = y1 + fraction * (y2 - y1)
        z = z1 + fraction * (z2 - z1)
        q = Quaternion([w,x,y,z])
        return q.unit_quat()

    def dot_product(self, q_param):
	"""

	Performs the dot product between two quaternion and returns the product in the form of a list of type [w,x,y,z].

	Arguments:
	    *q_param*
	        A list of type [w,x,y,z] used to represent a quaternion.

	"""

	w1, x1, y1, z1 = self.param
	w2, x2, y2, z2 = q_param
	return w1*w2 + x1*x2+ y1*y2 + z1*z2
	
    def to_rotation_matrix(self):
	"""

	Method to convert a quaternion to a rotation matrix.

	Return:
	    A rotation matrix
 
	"""

	a, b, c, d = self.param
	rot_mat = matrix([[a*a+b*b-c*c-d*d, 2*b*c-2*a*d, 2*a*c+2*b*d],
                          [2*b*c+2*a*d, a*a-b*b+c*c-d*d, 2*c*d-2*a*b],
                          [2*b*d-2*a*c, 2*a*b+2*c*d, a*a-b*b-c*c+d*d]])
	return rot_mat

    def check_for_mutations(self, mutationRate):
	"""

	Method to apply mutation operator on QuaternionGene based on mutationRate.

	Arguments:
	    *mutationRate*
	        Rate of mutation used to decide an application of the mutation operation.

	"""

        if uniform(0,1) < mutationRate:
            self.muate()
  
    def muate(self):
	"""

	Mutation operator to set a new value for self.param by randomly picking from the list of precomputed quaternion.

	"""

	range1 = self.ul_list[0]
	range2 = self.ul_list[1]

	ri_indx = randint(range1, range2)
	self.param = quat_vec[ri_indx]

    def __repr__(self):
	"""

	Returns the components of the quaternion objects as a string.

	"""

        return str(self.param)

class Genotype:
    """ 
    
    A class to store a collection of VectorGene and QuaternionGene that is used to represent the state of the components in the assembly. 
    
    """

    def __init__(self, gene_types):
	"""
	
	Constructor to initialise the Genotype.

	Arguments:
	    *gene_types*
	        A list of a list containing data to initialise the VectorGene and QuaternionGene which inturn initialise the Genotype object.

	"""

        self.fitness = 0
	self.nfitness = 0
	self.gof = 0
	self.clash = 0
	self.protruction = 0
	self.symres = 0
        self.genotype = []

       	for gene in gene_types:
	    t, ul_list, param = gene
            if t == 'f':
                self.genotype.append(FloatGene(pars[0], pars[1]))
            if t == 'i':
                self.genotype.append(IntGene(pars[0], pars[1]))
            if t == 'b':
                self.genotype.append(BinaryGene(pars[0]))
            if t == 's':
                self.genotype.append(SeqGene(pars[0]))
            if t == 'vt' or t == 'vr':
                self.genotype.append(VectorGene(t, ul_list,x,y,z))
            if t == 'cartesian':
                self.genotype.append(VectorGene(ul_list,param[0],param[1],param[2]))
            if t == 'quaternion':
                self.genotype.append(QuaternionGene(ul_list,param[0],param[1],param[2],param[3]))
    
    def get_fitness(self):
	"""
	
	Returns the fitness value of the genotype.

	"""

	return self.fitness

    def __getitem__(self, index):
	"""
	
	Returns the gene in the genotype based on the gene index.

	Arguments:
	    *index*
	        Index of the gene in the genotype.

	"""

        return self.genotype[index]

    def __len__(self):
	"""
	
	Returns the number of genes in the genotype.

	"""

        return len(self.genotype)

    #######def breed(self, no_of_iters, curr_iter, otherGenotype, mutationRate, crossRate):
    def breed(self, otherGenotype, mutationRate, crossRate):
	"""
	
	Generate a child genotype after applying crossover and mutation operation using two selected genotypes.

	Arguments:
	    *otherGenotype*
	        An instance of Genotype object.
	    *mutationRate*
	        Rate of mutation used.
	    *crossRate*
	        Rate of crossover used.

	Return:
	    Return a child genotype object.         

	"""

	ngene = self.__len__()
	child = self.uniform_crossover(otherGenotype, mutationRate, crossRate)	
	return child
	
    def uniform_crossover(self, otherGenotype, mutationRate, crossRate):
	"""
	
	Apply uniform crossover operation between two selected genotype.

	Arguments:
	    *otherGenotype*
	        An instance of Genotype object.
	    *mutationRate*
	        Rate of mutation used.
	    *crossRate*
	        Rate of crossover used.

	Return:
	    Return a child genotype object.         

	"""

        child = Genotype([])
	ngene = self.__len__()
	
	for ind, gene in enumerate(self.genotype):
		#print 'Gene inside multipoint crossover : ',gene
		newGene = gene.crossover(otherGenotype.genotype[ind], crossRate)
		newGene.check_for_mutations(mutationRate)
		child.genotype.append(newGene)
	return child
	
    def swap(self):
	"""
	
	Randomly select two genes (VectorGene and QuaternionGene) of two component in a genotype and swap them

	Return:
	    Return an instance of a genotype after swapping.        

	"""

	swap_geno = self.copy()
	
	n_subunit = self.__len__()/2
        ind = []
        for x in range(2):
            t = randrange(n_subunit)
            while t in ind:
                t = randrange(n_subunit)
            ind.append(t)
	unit1 = ind[0]+1
	unit2 = ind[1]+1
	unit_trans1 = 2*(unit1-1)
	unit_rot1 = 2*(unit1-1)+1
        unit_trans2 = 2*(unit2-1)
        unit_rot2 = 2*(unit2-1)+1
	#swap trans gene
	temp = swap_geno.genotype[unit_trans1]
	swap_geno.genotype[unit_trans1] = swap_geno.genotype[unit_trans2]
	swap_geno.genotype[unit_trans2] = temp

	#swap rot gene
        temp = swap_geno.genotype[unit_rot1]
        swap_geno.genotype[unit_rot1] = swap_geno.genotype[unit_rot2]
        swap_geno.genotype[unit_rot2] = temp

	return swap_geno

    def copy(self):
	"""
	
	Randomly select two genes (VectorGene and QuaternionGene) of two component in a genotype and swap them

	Return:
	    Return a copy of the genotype.         

	"""

        newGen = Genotype([])
        for x in self.genotype:
            newGen.genotype.append(x.copy())
        newGen.fitness = self.fitness
	#PAP addition
	newGen.gof = self.gof
	newGen.clash = self.clash
	newGen.protruction = self.protruction
	newGen.symres = self.symres
        return newGen

    def __repr__(self):
	"""
	
	Returns a string containing the details of the genotype including the gene content and its fitness values

	"""

        s = ''
        for x in self.genotype:
            s += str(x)+',\t'
        s += ': '+str(self.fitness)+', '+str(self.gof)+', '+str(self.clash)+', '+str(self.protruction)+', '+str(self.symres)
        return s

class Population:
    """ 
    
    A class to store a collection of Genotype objects. 
    
    """

    def __init__(self):
	"""
	
	Constructor to initialise the population object

	"""

        self.pop = []

    def __getitem__(self, index):
	"""
	
	Returns the genotype object in the population based on the genotype index.

	Arguments:
	    *index*
	        Index of the genotype in the population.

	"""

        return self.pop[index]

    def addGenotype(self, genotype):
	"""
	
	Method to append a genotype object to the population.

	"""

        self.pop.append(genotype)

    def size(self):
	"""
	
	Returns the number of genotypes in the population.

	"""

        return len(self.pop)

    def totalFitnessScore(self):
	"""
	
	Returns a string containing the total, average and the standard deviation of the fitness values in the population.

	"""

        total = 0
	avg = 0.0
	std = 0.0
	n = len(self.pop)
        for x in self.pop:
            total += x.fitness
	avg = total/n
	for x in self.pop:
	    std = std + (x.fitness - avg)**2
	std = sqrt(std/float(n-1))
        return str(total)+","+str(avg)+","+str(std)

    def pickSetOfBest(self, noOfBestGenotypes):
	"""
	
	Returns an instance of a population object containing n fittest genotypes.

	Arguments:
	    *noOfBestGenotype*
	        Number of fittest genotypes to pick.

	"""

        self.pop.sort(cmp = lambda x, y: cmp(-x.fitness, -y.fitness))
        newPop = []
        for g in range(noOfBestGenotypes):
            newPop.append(self.pop[g].copy())
        return newPop

    def pickSetOfWorst(self, noOfBestGenotypes):
	"""
	
	Returns an instance of a population object containing n worst fittest genotypes.

	Arguments:
	    *noOfBestGenotype*
	        Number of worst fittest genotypes to pick.

	"""

        self.pop.sort(cmp = lambda x, y: cmp(x.fitness, y.fitness))
        newPop = []
        for g in range(noOfBestGenotypes):
            newPop.append(self.pop[g].copy())
        return newPop

    def copy(self):
	"""
	
	Returns a copy of the population object.

	"""

	newPop = Population()
        for x in self.pop:
            newPop.addGenotype(x.copy())
	return newPop

    def pickBest(self):
	"""
	
	Returns the fittest genotype in the population.

	"""

        best = self.pop[0]
        for x in self.pop[1:]:
            if(x.fitness >= best.fitness):
                best = x
        return best

    def getBestScores(self):
	"""
	
	Returns a string containing the best fittest value in the pop and its corresponding components of the score (gof,clashscore)

	"""

        best = self.pop[0]
        for x in self.pop[1:]:
            if(x.fitness >= best.fitness):
                best = x
        scores = str(best.fitness)+","+str(best.gof)+","+str(best.clash)+","+str(best.protruction)+","+str(best.symres)
	return scores	

    def min_fitness(self):
	"""
	
	Returns the worst fittness score of the genotype in the population

	"""

        worst = self.pop[0]
        for x in self.pop[1:]:
            if(x.fitness <= worst.fitness):
                worst = x
	worstscore = worst.fitness
	return worstscore

    def max_fitness(self):
	"""
	
	Returns the best fittness score of the genotype in the population

	"""

        best = self.pop[0]
        for x in self.pop[1:]:
            if(x.fitness >= best.fitness):
                best = x
        bestscore = best.fitness
        return bestscore

    def avg_fitness(self):
	"""
	
	Returns the average fittness score of the genotypes in the population

	"""

	total = 0.0
        n = len(self.pop)
        for x in self.pop:
            total += x.fitness
        avg = total/n
	return avg

    def std_fitness(self):
	"""
	
	Returns the standard deviation of the fittness values of the genotypes in the population

	"""

	std = 0.0
	avg = self.avg_fitness()
        n = len(self.pop)
        for x in self.pop:
            std = std + (x.fitness - avg)**2
        std = sqrt(std/float(n-1))
	return std
 
    def pickWorst(self):
	"""
	
	Returns the worst fittness score of the genotype in the population

	"""

        worst = self.pop[0]
        for x in self.pop[1:]:
            if(x.fitness < worst.fitness):
                worst = x
        return worst
     
    #####def breed_1child(self, no_of_iters, curr_iter, mutation_const, crossRate, sel_method, sel_par):
    def breed_1child(self, mutation_const, crossRate, sel_method, sel_par):
	"""
	
	Returns a child genotype after applying crossover and mutation operation

	"""

        ma = sel_method(self, sel_par)
        pa = sel_method(self, sel_par)
	while ma == pa:
            pa = sel_method(self, sel_par)
	
	## Bread and return a single child
	#####child = ma.breed(no_of_iters, curr_iter, pa, mutation_const, crossRate)
	child = ma.breed(pa, mutation_const, crossRate)

	return child

    def breedNewPop(self, no_of_iters, curr_iter, mutation_const, crossRate, sel_method, sel_par):
	"""
	
	Returns a population object after performing a series of breeding opetation on the genotype in the current population

	"""

	newPop = Population()

	#Use it for general cases where components >= 2
	nbest = int(len(self.pop)*0.10)
	#For single component partial fit
        #nbest = 0
        
        for x in range(len(self.pop)-nbest):
            #####newPop.addGenotype(self.breed_1child(no_of_iters, curr_iter, mutation_const, crossRate, sel_method, sel_par))
            newPop.addGenotype(self.breed_1child(mutation_const, crossRate, sel_method, sel_par))
	
        top_best = self.pickBest()

        for i in range(nbest):
            #Get the randomly sawaped genotype and append it old_new_Pop        
	    tbest = deepcopy(top_best)
	    swap_geno = tbest.swap()
            top_mutated = Genotype([])
            for gene in swap_geno:
                newGene = gene.copy()
                newGene.check_for_mutations(0.1) # or mutation_const = 0.5
                top_mutated.genotype.append(newGene)
            newPop.addGenotype(top_mutated.copy())
	return newPop

    def __repr__(self):
	"""
	
	Returns a string containing the details of the population object

	"""

        rep = ''
        for x in self.pop:
            rep += str(x) + '\n'
        return rep

class Selection:
    """ 
    
    A class to implement the selection procedure used to select two genotype for the purpose of creating a child genotype. 
    
    """
    
    def roulette_wheel(self, pop):
	"""
	
	Return an instance of a Genotype object selected from the current population using Roulette wheel selection

	Arguments:
	    *pop*
	        Instance of a Population object.

	"""

	choice = uniform(0,1)
	currentChoice = 0
        for x in pop:
            currentChoice += x.nfitness
            if(currentChoice > choice):
                return x
    
    def tournament(self, pop, tourn_size=2):
	"""
	
	Return an instance of a Genotype object selected from the current population using Tournament selection

	Arguments:
	    *pop*
	        Instance of a Population object.
	    *tour_size*
	        Size of the tournament selection.

	"""

        tourn = []
        ind = []
        for x in range(tourn_size):
            t = randrange(pop.size())
            while t in ind:
                t = randrange(pop.size())
            ind.append(t)
            tourn.append(pop[t])
        tourn.sort(cmp = lambda x, y: cmp(-x.fitness, -y.fitness))
        return tourn[0]

    def SUS(self, pop):
        pass

    def local_select(self, pop):
        pass

    def migration_select(self, pop, no_of_subs):
        pass

#=====================================================================
# BELOW METHODS ARE NOT USED
#=====================================================================
    def breed_2child(self, no_of_iters, curr_iter, mutation_const, minmut, maxmut, sel_method, sel_par):
        ma = sel_method(self, sel_par)
        pa = sel_method(self, sel_par)
	while ma == pa:
            pa = sel_method(self, sel_par)

	## Bread and return two child
	child1, child2 = ma.breed(no_of_iters, curr_iter, pa, mutation_const, minmut, maxmut)
	return child1, child2

    def avg_nfitness(self):
	total = 0.0
        n = len(self.pop)
        for x in self.pop:
            total += x.nfitness
        avg = total/n
	return avg


    def std_nfitness(self):
	std = 0.0
	avg = self.avg_nfitness()
        n = len(self.pop)
        for x in self.pop:
            std = std + (x.nfitness - avg)**2
        std = sqrt(std/float(n-1))
	return std

    def normaliseScores(self):
        worstScore = self.pickWorst().fitness
        bestScore = self.pickBest().fitness
        if bestScore != worstScore:
            for x in self.pop:
                x.fitness -= worstScore
        else:
            for x in self.pop:
                x.fitness = 1

    def normalise_fitness(self):
	total,avg,std = self.totalFitnessScore()
	for x in self.pop:
	    x.fitness = x.fitness/total
	    
    def normalise_nfitness(self):
	total,avg,std = self.totalFitnessScore()
	min_fitness = self.min_fitness()
	max_fitness = self.max_fitness()
	for x in self.pop:
	    x.nfitness = (x.fitness-min_fitness)/(max_fitness-min_fitness)

    def rank_genotypes(self, sp=1.0):
	self.pop.sort(cmp = lambda x, y: cmp(-x.fitness, -y.fitness))
	for x in self.pop:
            rank = len(self.pop)-self.pop.index(x)
            x.fitness = 2-sp+2*(sp-1)*(rank-1.0)/(len(self.pop)-1)

    def breed_mutant(self, mutRate, minmut, maxmut):
	pop_size = len(self.pop)
	indx = sample(range(pop_size),1)[0]
	sel = self.pop[indx]
	mutant = sel.breed_mutant(mutRate, minmut, maxmut)
	return mutant

    def createMutant(self, mutRate, min_change, max_change):
	newPop = Population()
        for x in range(len(self.pop)):
            newPop.addGenotype(self.breed_mutant(mutRate, min_change, max_change))
	return newPop

#======================================================================

