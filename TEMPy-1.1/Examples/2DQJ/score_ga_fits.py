from TEMPy.Gamma import *
from cPickle import dump as dumpy
from TEMPy.ScoringFunctions import *
from itertools import *
#import matplotlib
#import pylab 
from numpy import *

def rotmat2axisang(mat):
	#print mat, mat[2,2]
	#exit(0)
	trace_mat = mat[0,0]+mat[1,1]+mat[2,2]
	theta1 = math.acos((trace_mat-1)/2)
	k = 1/(2*math.sin(theta1))
	n1 = k*(mat[2,1]-mat[1,2])
	n2 = k*(mat[0,2]-mat[2,0])
	n3 = k*(mat[1,0]-mat[0,1])
	return theta1,n1,n2,n3

def get_coor_mat(struct):
	coor_data = []
	n = len(struct.atomList)
	for atom in struct.atomList:
	    x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            m = atom.get_mass()
	    coor_data.append(x)
	    coor_data.append(y)
	    coor_data.append(z)
	    coor_data.append(m)
	xyzd = array(coor_data)
	xyzd = xyzd.reshape((n,4))
	xyzd = xyzd.T
	return xyzd

def get_coor_mat1(n_chain, selected_p_chain):
	ref_atoms, mov_atoms = get_equivalent_atoms(n_chain, selected_p_chain)	
        coor_data1 = []
	coor_data2 = []

	n = len(ref_atoms)
        for atom in ref_atoms:
            x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            m = atom.get_mass()
            coor_data1.append(x)
            coor_data1.append(y)
            coor_data1.append(z)
            coor_data1.append(m)
        xyzd1 = array(coor_data1)
        xyzd1 = xyzd1.reshape((n,4))
        xyzd1 = xyzd1.T

        for atom in mov_atoms:
            x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            m = atom.get_mass()
            coor_data2.append(x)
            coor_data2.append(y)
            coor_data2.append(z)
            coor_data2.append(m)
        xyzd2 = array(coor_data2)
        xyzd2 = xyzd2.reshape((n,4))
        xyzd2 = xyzd2.T
	return xyzd1, xyzd2
	
def get_transform(p_set1, p_set2):
	R = copy(p_set1[:, :])
	M = copy(p_set2[:, :])
	w = array([1.0], float)
	w = repeat(w, p_set1.shape[1])
	wR = array([R[0,:]*w, R[1,:]*w, R[2,:]*w])
	wM = array([M[0,:]*w, M[1,:]*w, M[2,:]*w])
	#print R.shape, M.shape
	mR = wR.mean(axis=1)
	mM = wM.mean(axis=1)
	
	#print mR, mM, mR.shape
	#print mR[0], mR[1], mR[2]
	#print wR, wM

	r1 = array(R[0,:]-mR[0])
	r2 = array(R[1,:]-mR[1])
	r3 = array(R[2,:]-mR[2])

	m1 = array(M[0,:]-mM[0])
	m2 = array(M[1,:]-mM[1])
	m3 = array(M[2,:]-mM[2])

	Rshifted = array([r1,r2,r3])
	Mshifted = array([m1,m2,m3])
	#print 'shape of rshift : ', Rshifted.shape
	#print Rshifted
	#exit(0)
	c = zeros((Rshifted.shape[0],Rshifted.shape[0]))	
	
	c[0,0] = dot(w,(Mshifted[0,:]*Rshifted[0,:]).T)
	c[0,1] = dot(w,(Mshifted[0,:]*Rshifted[1,:]).T)
	c[0,2] = dot(w,(Mshifted[0,:]*Rshifted[2,:]).T)

	c[1,0] = dot(w,(Mshifted[1,:]*Rshifted[0,:]).T)
	c[1,1] = dot(w,(Mshifted[1,:]*Rshifted[1,:]).T)
	c[1,2] = dot(w,(Mshifted[1,:]*Rshifted[2,:]).T)

	c[2,0] = dot(w,(Mshifted[2,:]*Rshifted[0,:]).T)
	c[2,1] = dot(w,(Mshifted[2,:]*Rshifted[1,:]).T)
	c[2,2] = dot(w,(Mshifted[2,:]*Rshifted[2,:]).T)


	#print c
	c = c / Rshifted.shape[1]
	u, s , v = linalg.svd(c)
	
	'''
	#Some checking on SVD calculations
	
	pinv = linalg.pinv(c)
	pin_svd = dot(dot(v.T,linalg.inv(diag(s))),u.T)

	print pinv
	print pin_svd
	'''

	R1 = dot(mat(v).H,u.T)

	if linalg.det(R1) < 0:
		print 'Det of rot mat is < 0'
		#print R1
		#print u
		#print s
		#print v
		exit(0)
		B = identity(3) 
		B[2,2] = linalg.det(mat(v).H*u.T)
		R1 = dot(dot(mat(v).H,B),u.T)
	
	#print 'Mean centered Ref and mov: ', mR, mM
	mM = mM.reshape((3,1))
	mR = mR.reshape((3,1))

	#print R1
	t = mR - dot(R1,mM)
	#print t, t.shape
	return R1,t


def superpose(p_set1, p_set2):
	r = 0
	t = 0
	maxIter = 50
	Rot = identity(3)
	TS = zeros((1,3))
	n = 0
	v_list = []
	for n in range(1):

 		R1, T1 = get_transform(p_set1, p_set2)

		T1 = T1.reshape((1,3))

		p_set2[:3,:] = dot(R1, p_set2[:3,:])

		p_set2[0,:] = p_set2[0,:] + T1[0,0]
		p_set2[1,:] = p_set2[1,:] + T1[0,1]
		p_set2[2,:] = p_set2[2,:] + T1[0,2]

		Rot = dot(R1,Rot)
		TS = dot(R1,TS.T).T
		TS[0,0] =  TS[0,0] + T1[0,0]
		TS[0,1] =  TS[0,1] + T1[0,1]
		TS[0,2] =  TS[0,2] + T1[0,2]
		
		#rot_points = p_set2[:,:]
		#v_list = get_vectors_1(rot_points)
		#vq_to_struct(v_list, output_file='After_rot'+str(n)+'.pdb')

		#print 'Current rot mat : ', R1
		#print 'New vec : ', TS
		#print 'New Rot : ', Rot
	return Rot, TS


def cal_ga_score(native_struct, predict_struct, emmap, resol, ncomp, template_grid, cvol, apix, asym_chains, sym_chains, rg_list):
	predict_assembly = Assembly(predict_struct.split_into_chains())
	predict_assembly.build_maps(resol,emmap)
	predict_assembly_map = predict_assembly.combine_maps()

	comps = predict_struct.split_into_chains()
	ncomp = len(comps)


	###predict_clash_score = round(scorer.get_sm_score(predict_struct, ncomp, template_grid, cvol, apix))
	#predict_protruction_penalty = scorer.get_protruction_score(predict, emmap, emmap.mean())
	###predict_protruction_penalty = 0.0
	###predict_mi = scorer.MI(predict_assembly_map,emmap)*ncomp
	#predict_mi = scorer.CCF(predict_assembly_map,emmap)*ncomp
	###predict_total_score = predict_mi + predict_clash_score + predict_protruction_penalty
	

	#comp_centroid = scorer.comp_centroid(native_struct, predict_struct, resolution)
	
	#print 'P_MI, P_clash, P_protruction, P_total, Centroid_score : ', predict_mi, predict_clash_score, predict_protruction_penalty, predict_total_score, comp_centroid

	#CALCULATE RMSD AND FRACTION_CORRECT SCORE

	whole_optimal_n_atom_list = []
	whole_optimal_p_atom_list = []
	total_chains = 0.0
	#cutoff for the centroid score
	cutoff = resolution/2
	num_good_comp = 0.0
	total_chain_rmsd = 0.0
	net_ang = 0.0
	net_trans = 0.0
	rmsd_list = []
	trans_list = []
	rot_list = []
	print 'Asymmetric chain calculation'
	for asym in asym_chains:
		total_chains = total_chains + 1.0
		print 'Chain ID : ',asym
		n_chain = native_struct.get_chain_ca(asym)
		p_chain = predict_struct.get_chain_ca(asym)
		rmsd = calculate_rmsd(n_chain, p_chain)
		total_chain_rmsd = total_chain_rmsd + rmsd
		#calculate centroid score
		#n_centroid = n_chain.calculate_centroid()
		#p_centroid = p_chain.calculate_centroid()
		n_centroid = n_chain.calculate_centre_of_mass()
		p_centroid = p_chain.calculate_centre_of_mass()

		dist = n_centroid.dist(p_centroid)
		
		#extreme = n_chain.get_extreme_values()
		#v1 = Vector(extreme[0],extreme[2],extreme[4])
		#v2 = Vector(extreme[1],extreme[3],extreme[5])
		#cutoff = v1.dist(v2)
		#cutoff = cutoff/2.0
		cutoff = rg_list[asym]
		if dist < cutoff:
			num_good_comp += 1.0
			print ('Component correctly placed based on Topology score (RMSD, native-predicted centroid distance, Rgyration cutoff): %6.2f %6.2f %6.2f' % (rmsd, dist, cutoff))
		else:
			print ('Component incorrectly placed based on Topology score (RMSD, native-predicted centroid distance, Rgyration cutoff): %6.2f %6.2f %6.2f' % (rmsd, dist, cutoff))
		

		rmsd_list.append(rmsd)
		#Add the chain atom list to build a whole structure for the optimal whole RMDS calculation
		for native_atom in n_chain: #uses __getitem__ function of the object
			whole_optimal_n_atom_list.append(native_atom)
		for predict_atom in p_chain: #uses __getitem__ function of the object
			whole_optimal_p_atom_list.append(predict_atom)


		#Call a function to calculate CPS between n_chain and selected_p_chain atoms
		#TO DO
		p_set1 = get_coor_mat(n_chain)

		com_ref = n_chain.calculate_centroid()
		com_mov = p_chain.calculate_centroid()

		#com_ref = n_chain.calculate_centre_of_mass()
		#com_mov = p_chain.calculate_centre_of_mass()


		offset = com_ref - com_mov
		p_chain.translate(offset.x,offset.y,offset.z)
		p_set2 = get_coor_mat(p_chain)

		r, t = superpose(p_set1, p_set2)
		ang, n1, n2, n3 = rotmat2axisang(r)
		
		print ('CPS score (Translation and rotation angle) : %6.2f, %6.2f' % (offset.mod(), math.degrees(ang)))
		net_ang = net_ang + math.degrees(ang)
		net_trans = net_trans + offset.mod()
		trans_list.append(offset.mod())
		rot_list.append(math.degrees(ang))
		p_chain.reset_position()


	sym_list = []
	for sym in sym_chains:
		#print sym
		sym_list.append(list(sym))

	#print sym_list

	print 'Symmetric chain calculation'
	for s in sym_list:
		#print 'Chain ID : ', s
		found_list = []
		for s1_chain_id in s:
			total_chains = total_chains + 1.0
			#print 'Chain ID : ', s1_chain_id
			#get the atom record (CA) of s_chain (predicted)
			n_chain = native_struct.get_chain_ca(s1_chain_id)
			min_rmsd = 10000
			for s2_chain_id in s:
				if s2_chain_id not in found_list:
					print s2_chain_id
					p_chain = predict_struct.get_chain_ca(s2_chain_id)
					#Collect the equivalent atoms (in case one chain has some missing atoms) and calculate the RMSD
					atom_list1, atom_list2 = get_equivalent_atoms(n_chain, p_chain)
					rmsd = calculate_rmsd(atom_list1, atom_list2)
					print rmsd
					#exit(0)
					if rmsd < min_rmsd:
						min_rmsd = rmsd
						chain_select = s2_chain_id
						selected_p_chain = p_chain
			total_chain_rmsd = total_chain_rmsd + min_rmsd	
			found_list.append(chain_select)
			print ('Chain pair and RMSD : %6.2f, %6.2f' % (s1_chain_id, chain_select, min_rmsd))
			rmsd_list.append(min_rmsd)
			atom_list1, atom_list2 = get_equivalent_atoms(n_chain, selected_p_chain)
			#Add the chain atom list to build a whole structure for the optimal whole RMDS calculation
			#for native_atom in n_chain:
                        for native_atom in atom_list1:
				whole_optimal_n_atom_list.append(native_atom)
			#for predict_atom in selected_p_chain:
                        for predict_atom in atom_list2:
				whole_optimal_p_atom_list.append(predict_atom)
		
			#calculate centroid score
			#n_centroid = n_chain.calculate_centroid()
			#p_centroid = selected_p_chain.calculate_centroid()

			n_centroid = n_chain.calculate_centre_of_mass()
			p_centroid = selected_p_chain.calculate_centre_of_mass()

			dist = n_centroid.dist(p_centroid)
			#print 'centroid dist : ', dist
 	                #extreme = n_chain.get_extreme_values()
        	        #v1 = Vector(extreme[0],extreme[1],extreme[2])
	           	#v2 = Vector(extreme[3],extreme[4],extreme[5])
	                #cutoff = v1.dist(v2)
			#cutoff = cutoff/2
 			cutoff = rg_list[s1_chain_id]
 			if dist < cutoff:
				num_good_comp += 1.0
				print ('Component correctly placed based on Topology score (RMSD, native-predicted centroid distance, Rgyration cutoff): %6.2f %6.2f %6.2f' % (rmsd, dist, cutoff))
			else:
				print ('Component incorrectly placed based on Topology score (RMSD, native-predicted centroid distance, Rgyration cutoff): %6.2f %6.2f %6.2f' % (rmsd, dist, cutoff))

				#print 'Number of components in place : ', num_good_comp

			#Call a function to calculate CPS between n_chain and selected_p_chain atoms
			#TO DO
			#p_set1 = get_coor_mat(n_chain)
			com_ref = n_chain.calculate_centroid()
			com_mov = selected_p_chain.calculate_centroid()
			#com_ref = n_chain.calculate_centre_of_mass()
			#com_mov = selected_p_chain.calculate_centre_of_mass()

			offset = com_ref - com_mov
			#print 'Translation vector : ', offset
			selected_p_chain.translate(offset.x,offset.y,offset.z)
			#p_set2 = get_coor_mat(selected_p_chain)
	
			p_set1, p_set2 = get_coor_mat1(n_chain, selected_p_chain)

			r, t = superpose(p_set1, p_set2)
			ang, n1, n2, n3 = rotmat2axisang(r)
			print ('CPS score (Translation and rotation angle) : %6.2f, %6.2f' % (offset.mod(), math.degrees(ang)))
			net_ang = net_ang + math.degrees(ang)
			net_trans = net_trans + offset.mod()
			trans_list.append(offset.mod())
			rot_list.append(math.degrees(ang))
			selected_p_chain.reset_position()

	avg_ang = net_ang/total_chains
	avg_trans = net_trans/total_chains

	#print 'Total chains : ', total_chains

	#Calculate the whole RMSD from the assembled chain atoms (with the information from the optimal symmetric chain paring)
	rmsd = calculate_rmsd(whole_optimal_n_atom_list, whole_optimal_p_atom_list)
	print ('Whole Ca RMSD : %6.2f' % rmsd)

	rmsd_avg = total_chain_rmsd/total_chains 
	#median_rmsd = median(rmsd_list)
	#median_trans = median(trans_list)
	#median_ang = median(rot_list)

	print ('Average chain rmsd = %6.2f' % rmsd_avg)
	#print 'Median rmsd = ', median_rmsd
	print ('Mean CPS_trans = %6.2f' % avg_trans)
	print ('Mean CPS_rot = %6.2f' % avg_ang)

	#print 'Median CPS_trans = ', median_trans
	#print 'Median CPS_rot = ', median_ang


	#Calculate centroid score
	fraction_good = num_good_comp/total_chains
	print ('Total number of correctly placed components : %3d' % int(num_good_comp))
	print ('Topology score : %6.2f' % fraction_good)


	#return predict_total_score, predict_mi, predict_clash_score, rmsd, rmsd_avg, median_rmsd, fraction_good, avg_ang, avg_trans, median_ang, median_trans
	#return rmsd, rmsd_avg, median_rmsd, fraction_good, avg_ang, avg_trans, median_ang, median_trans
	return rmsd, rmsd_avg, fraction_good, avg_ang, avg_trans

def get_equivalent_atoms(chain1, chain2):
	list1 = []
	list2 = []
	len_chain1 = len(chain1)
	len_chain2 = len(chain2)
	if len_chain1 < len_chain2:
		ref_chain = chain1
		compare_chain = chain2
	else:
		ref_chain = chain2
		compare_chain = chain1
	for atom1 in ref_chain:
		for atom2 in compare_chain:
			if atom1.res_no == atom2.res_no:
				list1.append(atom1)
				list2.append(atom2)
	#print 'LIST1 : ', list2
	#print 'LIST2 : ', list2
	#exit(0)		 
	return list1, list2

def calculate_rmsd(list1, list2):
	dists = []
	for a in range(len(list1)):
		#dists.append(list1[a].dist(list2[a])**2)
		dists.append(list1[a].distance_from_atom(list2[a])**2)
        dists = array(dists)

        #return dists.mean()
        return sqrt(dists.mean())


blurrer = StructureBlurrer()
scorer = ScoringFunctions()
# Read in EM map

emmap = MapParser.readMRC(sys.argv[1])
emmap.normalise()

native_struct = PDBParser.read_PDB_file('Test', sys.argv[2],hetatm=False,water=False)

predict_prefix = sys.argv[3]

native_assembly = Assembly(native_struct.split_into_chains())

resolution = float(sys.argv[4])

#Create a StructureBlurrer() instance - Used for clash score
blurrer = StructureBlurrer()

# apix used to record the atom position overlayed in the template map
apix = 3.5

# Get the volume occupied by the components in the assembly
cvol = emmap._get_component_volumes(native_struct, apix, blurrer) 

# Number of components in the assembly
comps = native_struct.split_into_chains()
ncomp = len(comps)

native_assembly.build_maps(resolution,emmap)
native_assembly_map = native_assembly.combine_maps()

#PAP ADDITION
template_grid = scorer.get_clash_map(emmap,apix)

#Read in the chain list for asymmetric subunits
asym_chains = list(sys.argv[5])
 
#Read in the chain list for symmetric subunits
sym_chains = sys.argv[6].split()

#get the individual chain structure and calculate the radius_of_gyration and store in a dictionary
chain_list = []
for c in asym_chains:
	chain_list.append(c)
	
for c1 in sym_chains:
	for c2 in c1:
		chain_list.append(c2)

#exit(0)
rg_list = {}
for c in chain_list:
	temp_chain = native_struct.get_chain(c)
	rg_list[c] = temp_chain.get_rgyration()
#print rg_list
#exit(0)

#native_struct.get_chain_ca(asym)

total_chains = 0

whole_optimal_n_atom_list = []
whole_optimal_p_atom_list = []

#cutoff for the centroid score
#cutoff = float(sys.argv[5])


#exit(0)

#Get the total score, mi and clash score for the set of predictions
total_score_list = []
mi_list = []
vo_list = []

gof = 0
f = open(predict_prefix+'.log','r')
with open(predict_prefix+'.log','r') as inF:
	for line in inF:
		if 'Number of GA Generations' in line:
			ngen = int(line.split(':')[1])
		if 'Number of GA runs' in line:
			nga = int(line.split(':')[1])
		if 'Mutual Information' in line:
			gof = 1
		if 'Cross Correlation Coefficient' in line:
			gof = 2

			
inF.close()

for i in range(nga):
	f = open(predict_prefix+'.log','r')
	with open(predict_prefix+'.log','r') as inF:
		for line in inF:
			if 'R'+str(i+1)+'G'+str(ngen) in line:
				temp = line.split(',')
				total_score_list.append(float(temp[1]))
				mi_list.append(float(temp[2]))
				vo_list.append(float(temp[3]))
				break
	inF.close()
	
#print total_score_list, mi_list, vo_list
#exit(0)

print 'Number of components : ',ncomp
print 'Component chain IDs : ', chain_list
print 'Asymmetric chain IDs : ',asym_chains
print 'Symmetric chain IDs : ',sym_chains


print '\nNATIVE SCORE'
print '--------------'
native_clash_score = round(scorer.get_sm_score(native_struct, ncomp, template_grid, cvol, apix))
#native_protruction_penalty = scorer.get_protruction_score(native, emmap, emmap.mean())
native_protruction_penalty = 0.0
if gof == 1:
	native_gof = scorer.MI(native_assembly_map,emmap)*ncomp
	print 'Goodness-of-fit score (Mutual Information) : ',native_gof
if gof == 2:
	native_gof = scorer.CCC(native_assembly_map,emmap)*ncomp
	print 'Goodness-of-fit score (Cross Correlation Coefficient) ',native_gof

native_total = native_gof + native_clash_score  + native_protruction_penalty
print 'Clash penalty : ',native_clash_score
print 'Total score : ',native_total


whole_rmsd_list = []
fraction_good_list = []
avg_rmsd_list = []
cps_rot = []
cps_trans = []
#median_rmsd = []
#median_ang = []
#median_trans = []
indx_list = range(1,nga+1)

print '\nDETAILED PREDICTION SCORES'
print '----------------------------'
print 'Total number of GA predicted fits in the directory ',predict_prefix,' : ', nga
#for i in range(1,no_predictions+1):
for i in range(1,nga+1):
	#i = 10
	print '\nPREDICTION SCORES FOR ',predict_prefix+'_'+str(i)+'.pdb'
	#for GA results
	#predict_struct = PDBParser.read_PDB_file(predict_prefix+'_'+str(i)+'.pdb')
	#for multifit results
        #predict_struct = PDBParser.read_PDB_file(predict_prefix+str(i)+'.pdb')
	predict_struct = PDBParser.read_PDB_file('Test',predict_prefix+'_'+str(i)+'.pdb', hetatm=False,water=False)
	#p_score_list.append(cal_ga_score(native_struct, predict_struct, emmap, resolution, ncomp, template_grid, cvol, apix, asym_chains, sym_chains))
	#rmsd, rmsd_avg, mid_rmsd, fraction_good, net_ang, net_trans, mid_ang, mid_trans = cal_ga_score(native_struct, predict_struct, emmap, resolution, ncomp, template_grid, cvol, apix, asym_chains, sym_chains, rg_list)
	rmsd, rmsd_avg, fraction_good, net_ang, net_trans = cal_ga_score(native_struct, predict_struct, emmap, resolution, ncomp, template_grid, cvol, apix, asym_chains, sym_chains, rg_list)

	#total_score_list.append(predict_total_score)
	#mi_list.append(mi)
	#vo_list.append(vo)
	whole_rmsd_list.append(rmsd)
	avg_rmsd_list.append(rmsd_avg)
	fraction_good_list.append(fraction_good)
	cps_rot.append(net_ang)
	cps_trans.append(net_trans)
	#median_rmsd.append(mid_rmsd)
	#median_ang.append(mid_ang)
	#median_trans.append(mid_trans)

	#exit(0)
#sorted_score_list = sorted(izip(total_score_list,mi_list,vo_list,whole_rmsd_list,avg_rmsd_list,median_rmsd,fraction_good_list,cps_rot,cps_trans,median_ang,median_trans,indx_list), reverse=True, key=lambda x: x[0])
#sorted_score_list = sorted(izip(total_score_list,mi_list,vo_list,whole_rmsd_list,avg_rmsd_list,fraction_good_list,cps_rot,cps_trans,indx_list), reverse=True, key=lambda x: x[0])
sorted_score_list = sorted(izip(indx_list,total_score_list,mi_list,vo_list,whole_rmsd_list,avg_rmsd_list,fraction_good_list,cps_rot,cps_trans), reverse=True, key=lambda x: x[1])


print '\nSUMMARY OF THE PREDICTIONS (RANKED IN ASSENDING ORDER)'
print '------------------------------------------------------'
print '\nGA_run index, Total score, goodness-of-fit*ncomp, Clash penalty, whole RMSD (A), Avg. RMSD (A), Topology score, CPS rotation (degrees), CPS translation (A)'
for i in range(0,nga):
	print ('%3d, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f ' % sorted_score_list[i])
		

