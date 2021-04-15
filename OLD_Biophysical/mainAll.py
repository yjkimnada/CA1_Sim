import numpy as np
import struct
import matplotlib.pyplot as plt
import time
import sys
import brian2 as br

import libcell as lb
import saveClass as sc

# this does not creates a new module, just executes the content of the file
# as if it was written in the same file!
execfile('sim_functs.py') 

#----------------------------------------------------------------------------
# Data saving object; # Data storage lists
data = sc.emptyObject()
data.vdata, data.vDdata, data.Gdata, data.Idata, data.stim = [], [], [], [], []

PAR1 = 1
PAR2 = 1

#----------------------------------------------------------------------------
# Simulation CONTROL
syntypes = ['alltree', 'local_regular', 'global_regular_lrand', 'global_regular_lreg']

data.model = 'CA1' # L23 or CA1
data.stimType = 'place' # poisson, place, replay, minis, nIter, SStim, DStim
data.synType = 'alltree' # synapse distribution - NULL, single, alltree, clustered, clust2, local_regular, global_regular_lrand, global_regular_lreg, clust2
# data.synType = syntypes[PAR1-1] # synapse distribution
data.actType = 'aDend' # passive, aSoma, aDend, active 
if (data.model == 'L23' and data.stimType == 'replay' and data.synType == 'clust2'):
	data.synType = 'clust3'

Nclusts = np.array([4])
Nsyn_per_clusts = np.array([60])

# 1min 20s / 10s simulation for active
# 30s / 10s simulation for passive
# 60s / 10s simulation for passive spines
data.modulateK = False # True or False, global modulation, similar to Losonczy 2006
data.modulateK_parents = False # True or False
data.modulateK_local = False # True or False
data.removeIspikes = 0 # 0-1. The probability of removing an inhibitory spike
data.modulateRmRa = False # True or False
data.modulateRmRaSeg = False # True or False
data.randomW = False # True or False
data.measureInputRes = False
data.iden = 5
data.randClust = False # randomize location of clustered synapses within the branch

data.constNMDA = False # if TRUE we use voltage independent NMDA receptors
data.modulateNa = True # True or False - switch off the Na in the branches

## only when data.stimType = 'place'
data.placeType = 'random_NR' # balanced, random_N (number of neurons), random_NR (number and peak firing rate)
if (data.placeType == 'balanced'): data.randomW = False

## only when data.stimType = 'nIter'
data.direction = 'IN' # 'OUT' or 'IN' only for nsyn - direction of stimulation sequence

data.AHS = False # active hotspots
## only when data.stimType = 'place'
data.stimseed = PAR2 # PAR2 - 1-10 

#---------------------------------------------------------------------------
data.SAVE = True
data.SHOWTRACES = False
data.SHOWSYNS = False

### number of iterations - only when data.stimType = 'place'
### this corresponds to different trials with identical synapses
data.nIter = 651 # max is usually 16
### time parameters
data.TSTOP = 10
if (data.model == 'L23'): data.TSTOP = 24
if (data.model == 'CA1' and data.stimType == 'replay'): data.TSTOP = 0.3

#----------------------------------------------------------------------------
# synapse parameters
#----------------------------------------------------------------------------
data.NMDA = True
data.GABA = True   
data.GABA_B = True

data.SPINES = True
data.Lmin = 60
weight_factor_A = 2 # multiplcative weight scale of the clustered synapses - AMPA
weight_factor_N = 2 # multiplcative weight scale of the clustered synapses - NMDA
data.g_factor = 1 # all synapses are scaled with this factor

## cluster parameters - only when data.stimType = 'place'
## clust: 240 = 1x240 = 4x60 = 12x20 = 48x5 = 240x1
## clust2: 240 = 4x60 = 8x30 = 12x20 = 24x10 = 48x5 = 120*2 = 240x1
Nclusts = np.array([4, 8, 12, 24, 48, 120, 240])
Nsyn_per_clusts = np.array([60, 30, 20, 10, 5, 2, 1])
#Nclusts = np.array([500, 1000, 2000]) # use this only if all synapses are clustered in CA1
#Nclusts = np.array([25]) # use this only if all synapses are clustered in CA1
# Nclusts = np.array([480, 960, 1920]) # use this only if all synapses are clustered in L23
#Nsyn_per_clusts = np.array([20])
# data.Lmin = 5

Nclust = Nclusts[PAR1-1]
Ncell_per_clust = Nsyn_per_clusts[PAR1-1]

if (data.stimType == 'minis'):
	Nclust = 10
	Ncell_per_clust = 1

mazeCenter = True
if (mazeCenter == True):
	inmazetype = 'midMaze'
else :
	inmazetype = 'randMaze'

execfile('init_params.py')

# #locDends = [1, 5, 7, 20, 21, 24, 32, 40, 48, 50, 52, 53, 56, 57, 59, 60, 64]
# locDends = [1, 5, 7, 20, 21, 32, 40, 48, 50, 52, 57]
# # locDends = [24, 53, 56, 59, 60, 64]
# data.locDend = [locDends[PAR1 - 1]]
# data.locDendRec = [data.locDend[0], 55] # for NIter - the stimulated branch
# data.x_DendRec = [ .8, .8]


#----------------------------------------------------------------------------
# Create neuron and add mechanisms
# if data.model == 'BS': model = lb.BS()


if data.model == 'L23': 
	model = lb.L23()
	if (data.synType == 'single'):
		model.dends[data.locDend[0]].nseg = 100

	if (data.modulateNa):
		model.gna_dend = 0

	if data.ACTIVE: 
		lb.init_active(model, axon=data.ACTIVE_soma,
							 soma=data.ACTIVE_soma, dend=data.ACTIVE_dend,
							 dendNa=data.ACTIVE_dend, dendCa=data.ACTIVE_dend)
		if data.ACTIVEhotSpot: lb.hotSpot(model)


if data.model == 'CA1': 
	model = lb.CA1()
	if (data.synType == 'single'):
		# model.dends[data.locDend[0]].nseg = 100
		if (data.modulateK == True): # global modulation - similar to Losonczy et al. 2006
			model.gka = model.gkdr / 2     
			model.gka_trunk = model.gka_trunk / 2    

		if (data.modulateRmRa == True):
		  model.dends[data.locDend[0]].g_pas = model.dends[data.locDend[0]].g_pas / 4 
		  model.dends[data.locDend[0]].Ra = model.dends[data.locDend[0]].Ra / 2   

		if (data.modulateRmRaSeg == True):
		  sec = model.dends[data.locDend[0]]
		  i = 0
		  for seg in sec:
		      i = i + 1
		      if ((i > (data.locSeg[0]*100 - 5)) & (i < (data.locSeg[1] * 100+5))):
		          seg.g_pas = seg.g_pas / 4 / np.sqrt(2)
		          seg.diam = seg.diam * np.sqrt(2)    

	if (data.modulateNa):
		model.gna = 0     


	if data.ACTIVE: 
		lb.init_activeCA1(model, soma=data.ACTIVE_soma, dend=data.ACTIVE_dend)
		if data.ACTIVEhotSpot: 
			for d in data.hotspot_branches:
				nseg_d = model.dends[d].nseg
				model.dends[d].nseg = max(nseg_d, 10)
				lb.CA1_hotSpot(model, d, 10)
	else :
		lb.init_passiveCA1(model)


if data.model == 'CA1': 
	if (data.synType == 'single'):
		if (data.modulateK_local == True):
			gfactor = 0
			model.dends[data.locDend[0]].gkabar_kap = gfactor * model.dends[data.locDend[0]].gkabar_kap
			model.dends[data.locDend[0]].gkabar_kad = gfactor * model.dends[data.locDend[0]].gkabar_kad
			model.dends[data.locDend[0]].gkdrbar_kdr = gfactor * model.dends[data.locDend[0]].gkdrbar_kdr 

		if (data.modulateK_parents == True):
			for d in data.parents:
				gfactor = 0
				model.dends[d].gkabar_kap = gfactor * model.dends[d].gkabar_kap
				model.dends[d].gkabar_kad = gfactor * model.dends[d].gkabar_kad
				model.dends[d].gkdrbar_kdr = gfactor * model.dends[d].gkdrbar_kdr 



if (data.ICLAMP):
	if data.iclampLoc[0]=='soma':
		lb.add_somaStim(model, data.iclampLoc[1], onset=data.iclampOnset,
						dur=data.iclampDur, amp=data.iclampAmp)
	if data.iclampLoc[0]=='dend':
		lb.add_dendStim(model, data.iclampLoc[1], data.iclampLoc[2],
				 onset=data.iclampOnset, dur=data.iclampDur, amp=data.iclampAmp)

if (data.measureInputRes == True):
	lb.add_dendStim4(model, dends=data.locDendRec, onset=data.iclampOnset, dur=data.iclampDur, amp=data.iclampAmp)





#----------------------------------------------------------------------------
# Generate synapse locations
#----------------------------------------------------------------------------
# nsyn_soma = int(data.Insyn / 2)
nsyn_soma = 80 # 80
if (data.stimType == 'minis'): nsyn_soma = 0
trunk_ids = np.array([-1, 76, 80, 88, 90, 94, 102, 104, 106, 110])
# trunk_lengths = np.array([60, 57, 21, 7,  31, 11, 33, 7, 5, 20]) # total 192 um
trunk_lengths = np.array([60, 45, 16, 5,  21, 8, 22, 4, 3, 10]) # decreasing density towards the tuft
trunk_probs = trunk_lengths / float(sum(trunk_lengths))
trunk_ids2 = np.array([76, 80, 88, 90, 94, 102, 104, 106, 110, 112, 114, 118]) # for measuring EPSP properties of synapses along the trunk, 53-307 um

data.ind_clust = np.arange(0, data.Ensyn)
recNMDA = [880, 881, 882, 883]


np.random.seed(data.stimseed)

if (data.synType != 'NULL'):
	print (data.Ensyn)
	# Elocs: synapse locations - list, with elements [dendrite number, synapse location]
	if (data.synType == 'alltree'):
		data.Elocs = genRandomLocs(data.Ensyn, seed=100*data.stimseed+1)
		np.random.shuffle(data.Elocs)
		print ('number of RANDOM excitatory synapses:', data.Ensyn)
	elif (data.synType == 'local_regular'):
		data.Elocs = genLocalRegularLocs(data.Ensyn, seed=100*data.stimseed+1, regular_location=True)
		print ('number of RANDOM excitatory synapses:', data.Ensyn)
	elif (data.synType == 'trunk'):
		data.Elocs = genRandomLocs(data.Ensyn, seed=100*data.stimseed+1, dend_ids=trunk_ids2)           
	# data.Elocs = genRandomLocs(data.Ensyn, 25) # number of synapses, 
		print ('number of RANDOM excitatory synapses on the TRUNK:', data.Ensyn)

	elif (data.synType == 'global_regular_lrand'):
		# data.Elocs = mainAll.genAllLocs(5.28, [0, 1]) # distance between synapses - 2000 synapse generated
		if data.model == 'L23': 
			data.Elocs = genUnifLocsL23(2.92, 1920, localRandom=True) # distance between synapses - 2000 synapse generated
		else:
			data.Elocs = genUnifLocs(4.877, 2000, Ngroups=40, prime=9, localRandom=True) # distance between synapses - 2000 synapse generated
		data.Ensyn = len(data.Elocs)
		print ('number of globally UNIFORM, locally RANDOM excitatory synapses:', data.Ensyn)

	elif (data.synType == 'global_regular_lreg'):
		# data.Elocs = mainAll.genAllLocs(5.28, [0, 1]) # distance between synapses - 2000 synapse generated
		if data.model == 'L23': 
			data.Elocs = genUnifLocsL23(2.92, 1920) # distance between synapses - 2000 synapse generated
		else:
			data.Elocs = genUnifLocs(4.877, 2000, Ngroups=40, prime=9) # distance between synapses - 2000 synapse generated
		data.Ensyn = len(data.Elocs)
		print ('number of globally UNIFORM, locally REGULAR excitatory synapses:', data.Ensyn)

	elif (data.synType == 'somatic'):
		data.Elocs = genLocalRegularLocs(data.Ensyn, seed=100*data.stimseed+1)
		for i in range(data.Ensyn):
			data.Elocs[i][0] = 50
			data.Elocs[i][1] = 0.5
		print ('number of somatic excitatory synapses:', data.Ensyn)
	elif (data.synType == 'clustered'):
		data.Elocs, data.ind_clust = genClustLocs(data.Ensyn, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter) # number of synapses, 
		print ('number of CLUSTERED excitatory synapses:', len(data.Elocs))
		print (Nclust, 'cluster of size', Ncell_per_clust, 'generated')
	elif (data.synType == 'clust2'):
		data.Elocs, data.ind_clust = addClustLocs(data.Ensyn, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter, clocs=data.cdends[PAR2 - 1], Lmin=data.Lmin) # number of synapses,
		recNMDA = [880 + Ncell_per_clust / 2, 880 + Ncell_per_clust + Ncell_per_clust / 2, 880 + Ncell_per_clust * 2 + Ncell_per_clust / 2, 880 + Ncell_per_clust * 3 + Ncell_per_clust / 2]
		print ('number of CLUSTERED excitatory synapses:', len(data.Elocs))
		print (Nclust, 'cluster of size', Ncell_per_clust, 'generated')
	elif (data.synType == 'clust3'):
		# Elocs, data.ind_clust = addClustLocs(data.Ensyn, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter) # number of synapses,
		Elocs, data.ind_clust = addClustLocs(data.Ensyn, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter, clocs=data.cdends[PAR2 - 1]) # number of synapses,
		DSyn = lb.synDist(model, Elocs)
		data.Elocs = rearrangeBgLocs(Elocs, DSyn, data.ind_clust, range(data.Ensyn-Nclust*Ncell_per_clust,data.Ensyn))

		# Elocs, data.ind_clust = addClustLocs(data.Ensyn-Nclust*Ncell_per_clust, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter) # number of synapses,
		# data.Elocs = addBgLocs(Elocs, data.ind_clust)
		# Elocs, data.ind_clust = addClustLocs(2000, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=mazeCenter) # number of synapses,
		# data.Elocs = addBgLocs(Elocs, range(880,1120))
		# data.Elocs = Elocs
		recNMDA = [880 + Ncell_per_clust / 2, 880 + Ncell_per_clust + Ncell_per_clust / 2, 880 + Ncell_per_clust * 2 + Ncell_per_clust / 2, 880 + Ncell_per_clust * 3 + Ncell_per_clust / 2]
		print ('number of CLUSTERED excitatory synapses:', len(data.Elocs))
		print (Nclust, 'cluster of size', Ncell_per_clust, 'generated')

	elif (data.synType == 'single'):
		data.Elocs = genDendLocs(data.locDend, data.Ensyn, data.locSeg)

	# data.locDendRec = [8, 59, 64, 71]#, 13, 33, 108] 
	# data.xDendRec = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
	# data.locDendRec = [data.Elocs[950][0], data.Elocs[1050][0], 52, 28] 
	# data.locSpineRec = [950, 1050, 75, 134] 

	# data.locDendRec = [data.Elocs[569][0], data.Elocs[77][0], data.Elocs[950][0], data.Elocs[1050][0]] # pre spikes at 4400-5500
	# data.xDendRec = [data.Elocs[569][1], data.Elocs[77][1], data.Elocs[950][1], data.Elocs[1050][1]] 
	# data.locSpineRec = [569, 77, 950, 1050]

	# for j_syn in np.arange(0, 240):
	#     data.Elocs[j_syn+880] = [92, j_syn / 240.]
		# data.Elocs[j_syn+880] = [92, 0.1]
		# data.Elocs[j_syn+880] = [92, 0.9]
	# 
	# data.locDendRec = data.cdends[PAR2 - 1] #, 13, 33, 108] 
	# data.xDendRec = [0.5, 0.5, 0.5, 0.5]#, 0.5, 0.5, 0.5, 0.5] 
	data.locDendRec = [5, 59, 6, 6] #, 13, 33, 108] 
	data.xDendRec = [0.5, 0.6, 0.89, 0.75]#, 0.5, 0.5, 0.5, 0.5] 

	if (data.randClust):
		for i in range(240):
			data.Elocs[i+880][1] = np.random.uniform()


	print(len(data.Elocs))

	if data.GABA:
		Isomalocs = []
		np.random.seed(data.stimseed + 10000)
		for p in np.arange(0,nsyn_soma):
			if (data.model == 'CA1'): # L23 or CA1
				id_den = np.random.choice(trunk_ids, 1, p=trunk_probs)[0]
				loc_syn = np.random.rand(1)[0]
				Isomalocs.append([id_den, loc_syn])
				# Isomalocs.append([-1, 0.5])
			else:
			   Isomalocs.append([-1, 0.5])# -1 is not an index - it means that the synapse is at the soma
		Idendlocs = genRandomLocs(int(data.Insyn - nsyn_soma), data.stimseed + 10000)
		np.random.shuffle(Idendlocs)
		np.random.shuffle(Isomalocs)
		data.Ilocs = Isomalocs + Idendlocs

		if (data.synType == 'somatic'):
			for i in range(data.Insyn):
				data.Ilocs[i][0] = 50
				data.Ilocs[i][1] = 0.5

		print ('number of inhibitory synapses:', len(data.Ilocs))


	# for isyn in np.arange(len(data.Elocs)):
	#   data.Elocs[isyn] = [0, 0.5]
	# for isyn in np.arange(len(data.Ilocs)):
	#   data.Ilocs[isyn] = [-1, 0.5]

	lb.add_syns(model, data)

	if ((data.synType == 'clustered') + (data.synType == 'clust2') + (data.synType == 'clust3')):
		if (weight_factor_N != 1):
			for i in data.ind_clust:
				model.ncNMDAlist[i].weight[0] = weight_factor_N * model.ncNMDAlist[i].weight[0] 
		if (weight_factor_A != 1):
			for i in data.ind_clust:
				model.ncAMPAlist[i].weight[0] = weight_factor_A * model.ncAMPAlist[i].weight[0] 

	if (data.randomW == True):
		w = model.ncAMPAlist[0].weight[0]
		wN = model.ncNMDAlist[0].weight[0]
		print('old weight', w)
		for i in np.arange(data.Ensyn):
			model.ncAMPAlist[i].weight[0] = np.random.uniform(w*0.25, w*1.75)
			model.ncNMDAlist[i].weight[0] = np.random.uniform(wN*0.25, wN*1.75)
			if (i < 10): 
				print('new weight:', model.ncAMPAlist[i].weight[0])
	
	else:
		print('synaptic weights are CONSTANT')


# for i in np.arange(data.Ensyn):
# 	model.ncAMPAlist[i].weight[0] = 0
# 	model.ncNMDAlist[i].weight[0] = 0

# model.ncAMPAlist[914].weight[0] = data.Agmax/1000.
# model.ncNMDAlist[914].weight[0] = data.Ngmax/1000.

# for i in np.arange(data.Insyn):
# 	model.ncGABAlist[i].weight[0] = 0
# 	model.ncGABA_Blist[i].weight[0] = 0

# #----------------------------------------------------------------------------
# run simulation
#----------------------------------------------------------------------------
print ('synapse distribution: ', data.synType)
print ('stimulus type:', data.stimType)
print ('cell type:', data.actType)
print ('simulation length:', data.TSTOP, 's')
print ('gAMPA:', data.Agmax)


if (data.stimType == "poisson") :
# def SIM_PoissonIteration(Ensyn, Erate,  Insyn, Irate, duration, nIter):
	SIM_PoissonIteration(data.Ensyn, data.Erate, data.Insyn, data.Irate, data.TSTOP, data.nIter)
elif (data.stimType == "place"):
# def SIM_PlaceIteration(type, Ensyn, Insyn, Erate, Irate, duration, Nrep, stimseed):
	SIM_PlaceIteration(data.placeType, data.Ensyn, data.Insyn, data.Erate, data.Irate, data.TSTOP, data.nIter, data.stimseed, recNMDA, data.removeIspikes)
elif (data.stimType == "replay"):
# def SIM_PlaceIteration(type, Ensyn, Insyn, Erate, Irate, duration, Nrep, stimseed):
	SIM_ReplayIteration(data.placeType, data.Ensyn, data.Insyn, data.Erate, data.Irate, data.TSTOP, data.nIter, recNMDA)
elif (data.stimType == "minis"):
# def SIM_minisDistribution(data):
	SIM_minisDistribution(data)
elif (data.stimType == "nIter"):
# def SIM_nsynIteration(maxNsyn, tInterval, onset, direction='OUT'):
	SIM_nsynIteration(data.Ensyn, data.tInterval, data.st_onset, direction=data.direction)
elif ((data.stimType == "SStim") + (data.stimType == "DStim")):
# def SIM_currentSteps(iRange):
	SIM_currentSteps(data.iRange)

# # ----------------------------------------------------------------------------
# # show the traces

if data.SHOWTRACES:
	print ('plotting traces ...')
	import cell_traces as ct
	if ((data.stimType == 'minis') + (data.stimType == 'nIter') + (data.stimType=='SStim')):
		# ct.plotTraces(data, data.TSTOP * 1000, len(data.locDendRec))
		# ct.plotV_G(data, model.seclength)
		ct.plotV_dV(data)
	elif (data.stimType=='DStim'):
		ct.plotTraces(data, data.TSTOP * 1000, len(data.locDendRec))
		ct.plotAPprop(data)
		ct.plotV_dV(data)
		# ct.plotV_G(data, model.seclength)
	else :
		ct.plotResp(data, data.TSTOP * 1000)


# #----------------------------------------------------------------------------
# # Save data - number of clusters are encoded in the locBias parameter

if data.SAVE:
	print ('saving data ...')
	if ((data.synType == 'clust2') + (data.synType == 'clust3')):
		if (data.SPINES):
			outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '_' + inmazetype + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'spines'
			if (data.ACTIVEhotSpot):
				outdir = outdir + '_hotspot'
			if (data.randClust):
				outdir = outdir + '_randClust'
		else :
			outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '_' + inmazetype + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'nospines'                   
			if (data.removeIspikes > 0):
				outdir = outdir + '_Ireduced' + str(data.removeIspikes)
			if (data.ACTIVEhotSpot):
				outdir = outdir + '_hotspot'

	else :
		outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType
		if (data.constNMDA == True):
			outdir = outdir + '_NMDAc'

	if (data.stimType == 'nIter'):
		outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '/dend' + str(data.locDend[0])
	if ((data.stimType == 'SStim') + (data.stimType == 'DStim')):
		outdir = './' + data.model + '/' + data.stimType + '/' + data.actType
	import cell_save as cs
	cs.save_sim(data, out_binary=True, out_vdend=True, out_pickle=False, outdir=outdir, dt_save=0.2)

# #----------------------------------------------------------------------------
# # show the synapses

if data.SHOWSYNS:
	print ('plotting the cell with the synapses ...')
	import cell_draw as cd
	if ((data.stimType=='nIter')) :
		cd.plot_syns(data, model, False)
	else : 
		if (data.stimType == 'replay'): 
			iB = range(data.Ensyn-Nclust*Ncell_per_clust,data.Ensyn)
		else :
			iB = []
		cd.plot_syns(data, model, True, iB)

# #----------------------------------------------------------------------------
# # inter-synapse matrix

# # print 'calculating the inter-synapse distance matrix'

# DSyn = lb.synDist(model, data.Elocs)
# # imgplot = plt.imshow(DSyn)
# # plt.show(block=False)


# # mainAll.lb.h.shead[12].g_pas
# # model.dends[12].g_pas

#########################################################
#########################################################
## saving dendritic responses

# data = mainAll.data
# model = mainAll.model
# inmazetype = mainAll.inmazetype
# Nclust = mainAll.Nclust
# Ncell_per_clust = mainAll.Ncell_per_clust
# weight_factor_A = mainAll.weight_factor_A
# weight_factor_N = mainAll.weight_factor_N

# import numpy as np
# import cell_save as cs


# outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/D' + data.actType + '_' + inmazetype + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'spines'

# if (data.stimType == 'place'):
#   print('saving average vdend - place')
#   bfname = './'+outdir+'/mean_VD_I_'+str(data.stimseed)+'.bin'
#   print (bfname)
#   cs.save_ave_place(data.vDdata, data.nIter, bfname)

#   bfname = './'+outdir+'/mean_GN_'+str(data.stimseed)+'.bin'
#   cs.save_ave_place(data.Gdata, data.nIter, bfname)

#   bfname = './'+outdir+'/mean_IN_'+str(data.stimseed)+'.bin'
#   cs.save_ave_place(data.Idata, data.nIter, bfname)

# if (data.stimType == 'replay'):
#   print('saving average vdend - replay')
#   bfname = outdir+'/mean_VD_I'+str(data.stimseed)+'.bin'
#   print (bfname)
#   cs.save_ave_replay(data.vDdata, data.nIter, 20, bfname)

#   bfname = outdir+'/mean_GN_'+str(data.stimseed)+'.bin'
#   cs.save_ave_replay(data.Gdata, data.nIter, 20, bfname)

#   bfname = outdir+'/mean_IN_'+str(data.stimseed)+'.bin'
#   cs.save_ave_replay(data.Idata, data.nIter, 20, bfname)


############################################################
# brs = np.arange(151) - 100
# counts = np.zeros((150, 80), dtype=int)
# # counts = np.zeros((150, 12), dtype=int)

# for i_trial in range(data.nIter):
#     vv = data.vDdata[i_trial]
#     k = 0
#     for i_dendrite in range(4):
#         vvv = vv[i_dendrite]
#         mv = np.reshape(vvv[0:50000], (20, 2500))
#         # mv = np.reshape(vvv, (20, 1501))
#         for i_startpoint in range(20):
#             # hvd = np.histogram(mv[i_startpoint,550:1000], brs)
#             hvd = np.histogram(mv[i_startpoint,:], brs)
#             counts[:,k] = counts[:,k] + hvd[0]
#             k = k + 1

#         # mv = np.reshape(vvv, (3, 1501))
#         # for i_startpoint in range(3):
#         #     hvd = np.histogram(mv[i_startpoint,550:1000], brs)
#         #     counts[:,k] = counts[:,k] + hvd[0]
#         #     k = k + 1

# outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '_' + inmazetype + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'nospines'
# bfname = './'+outdir+'/hist_vD4_'+str(data.stimseed)+'.bin'
# binfile = file(bfname, 'wb')
# # and write out two integers with the row and column dimension
# header = struct.pack('2I', counts.shape[0], counts.shape[1])
# binfile.write(header)
# # then loop over columns and write each
# for i in range(counts.shape[1]):
#     ddata = struct.pack('%iI' % counts.shape[0], *counts[:,i])
#     binfile.write(ddata)
# binfile.close()

# fromSection = model.soma
# mainAll.lb.h.distance(0, 0.5, sec=fromSection)
# toSection = model.dends[57]
# mainAll.lb.h.distance(0.5, sec=toSection)
