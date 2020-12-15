import numpy as np
import struct
import matplotlib.pyplot as plt
import time
import sys
# import brian as br

import libcell as lb
import saveClass as sc


PAR1 = 1 # controls the level of clustering
PAR2 = 0 # controls the random seed

#----------------------------------------------------------------------------
# Data saving object; # Data storage lists
data = sc.emptyObject()
data.vdata, data.vDdata, data.Gdata, data.Idata, data.stim = [], [], [], [], []

#----------------------------------------------------------------------------
# Simulation CONTROL
data.SAVE = True # save the results of the simulations into file
data.SHOWTRACES = False # plot the input and the dendritic and somatic voltage traces
data.SHOWSYNS = False # show the morphology of the cell with the location of the synapses

### number of iterations - only when data.stimType = 'place'
### this corresponds to different trials with identical synapses
data.nIter = 16 # max is usually 16, may not work propely with < 2
data.TSTOP = 10

#---------------------------------------------------------------------------
# biophysics of dendrites and soma
data.actType = 'aDend' # passive, aSoma, aDend, active 
data.modulateNa = False # True or False - switch off the Na in the branches, otherwise keep them constant
data.AHS = False # active hotspots - when you have 12 clusters, you can add them on branches with Na hotspots

data.stimseed = PAR2+1 # PAR2 - we only have stimulus for PAR2=0 :-(

#---------------------------------------------------------------------------
# connectivity parameters
Nclusts = np.array([4, 8, 12, 24, 48, 120, 240])
Nsyn_per_clusts = np.array([60, 30, 20, 10, 5, 2, 1])

Nclust = Nclusts[PAR1]
Ncell_per_clust = Nsyn_per_clusts[PAR1]


# execfile('init_params.py')
exec(open("./init_params.py").read())
exec(open("./sim_functs.py").read())

#----------------------------------------------------------------------------
# Create neuron and add mechanisms
# if data.model == 'BS': model = lb.BS()

model = lb.CA1()

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

#----------------------------------------------------------------------------
# Generate synapse locations
#----------------------------------------------------------------------------
# nsyn_soma = int(data.Insyn / 2)
nsyn_soma = 80 # 80
trunk_ids = np.array([-1, 76, 80, 88, 90, 94, 102, 104, 106, 110])
trunk_lengths = np.array([60, 45, 16, 5,  21, 8, 22, 4, 3, 10]) # decreasing density towards the tuft
trunk_probs = trunk_lengths / float(sum(trunk_lengths))

data.ind_clust = np.arange(0, data.Ensyn)


np.random.seed(data.stimseed)

if (data.synType != 'NULL'):
    print (data.Ensyn)
    if (data.synType == 'clust2'):
        data.Elocs, data.ind_clust = addClustLocs(data.Ensyn, Nclust, Ncell_per_clust, 25 + data.stimseed, midle=True, clocs=data.cdends[PAR2], Lmin=data.Lmin) # number of synapses,
        print ('number of CLUSTERED excitatory synapses:', len(data.Elocs))
        print (Nclust, 'cluster of size', Ncell_per_clust, 'generated')
    elif (data.synType == 'single'):
        data.Elocs = genDendLocs(data.locDend, data.Ensyn, data.locSeg)

    print(len(data.Elocs))

    if data.GABA:
        Isomalocs = []
        np.random.seed(data.stimseed + 10000)
        for p in np.arange(0,nsyn_soma):
            # id_den = np.random.choice(trunk_ids, 1, p=trunk_probs)[0]
            # loc_syn = np.random.rand(1)[0]
            # Isomalocs.append([id_den, loc_syn])
            Isomalocs.append([-1, 0.5])
        Idendlocs = genRandomLocs(int(data.Insyn - nsyn_soma), data.stimseed + 10000)
        np.random.shuffle(Idendlocs)
        np.random.shuffle(Isomalocs)
        data.Ilocs = Isomalocs + Idendlocs


        print ('number of inhibitory synapses:', len(data.Ilocs))

    lb.add_syns(model, data)

    if (data.synType == 'clust2'):
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

"""
synfilename = './CA1stims/synlocs_clust4_rep16_stimseed1.bin'
f=open(synfilename, "rb")
print(struct.unpack('2I', f.read(8)))
for ii in np.arange(2000):
    data.Elocs[ii][0] = int(struct.unpack('d', f.read(8))[0])
for ii in np.arange(200):
    data.Ilocs[ii][0] = int(struct.unpack('d', f.read(8))[0])
for ii in np.arange(2000):
    data.Elocs[ii][1] = struct.unpack('d', f.read(8))[0]
for ii in np.arange(200):
    data.Ilocs[ii][1] = struct.unpack('d', f.read(8))[0] - 1
f.close()
"""

# for i in np.arange(data.Ensyn):
#   model.ncAMPAlist[i].weight[0] = 0
#   model.ncNMDAlist[i].weight[0] = 0

# model.ncAMPAlist[914].weight[0] = data.Agmax/1000.
# model.ncNMDAlist[914].weight[0] = data.Ngmax/1000.

# for i in np.arange(data.Insyn):
#   model.ncGABAlist[i].weight[0] = 0
#   model.ncGABA_Blist[i].weight[0] = 0

# #----------------------------------------------------------------------------
# run simulation
#----------------------------------------------------------------------------
print ('synapse distribution: ', data.synType)
print ('stimulus type:', data.stimType)
print ('cell type:', data.actType)
print ('simulation length:', data.TSTOP, 's')
print ('gAMPA:', data.Agmax)


if (data.stimType == "place"):
# def SIM_PlaceIteration(type, Ensyn, Insyn, Erate, Irate, duration, Nrep, stimseed):
    SIM_PlaceIteration(data.placeType, data.Ensyn, data.Insyn, data.Erate, data.Irate, data.TSTOP, data.nIter, data.stimseed, data.removeIspikes)
elif (data.stimType == "minis"):
# def SIM_minisDistribution(data):
    SIM_minisDistribution(data)
elif ((data.stimType == "SStim") + (data.stimType == "DStim")):
# def SIM_currentSteps(iRange):
    SIM_currentSteps(data.iRange)

# # ----------------------------------------------------------------------------
# # show the traces

if data.SHOWTRACES:
    print ('plotting traces ...')
    import cell_traces as ct
    if (data.stimType == 'minis'):
        # ct.plotTraces(data, data.TSTOP * 1000, len(data.locDendRec))
        # ct.plotV_G(data, model.seclength)
        ct.plotV_dV(data)
    else :
        ct.plotResp(data, data.TSTOP * 1000)


# #----------------------------------------------------------------------------
# # Save data

if data.SAVE:
    print ('saving data ...')
    if (data.synType == 'clust2'):
        if (data.SPINES):
            outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'spines'
            if (data.ACTIVEhotSpot):
                outdir = outdir + '_hotspot'
        else :
            outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType + '_' + str(Nclust) + 'clusterOf' + str(Ncell_per_clust) + '_wA' + str(weight_factor_A) + '_wN' + str(weight_factor_N) + 'nospines'                   
            if (data.removeIspikes > 0):
                outdir = outdir + '_Ireduced' + str(data.removeIspikes)
            if (data.ACTIVEhotSpot):
                outdir = outdir + '_hotspot'

    else :
        outdir = './' + data.model + '/' + data.synType + '/' + data.stimType + '/' + data.actType

    import cell_save as cs
    cs.save_sim(data, out_binary=True, out_vdend=True, out_pickle=False, outdir=outdir, dt_save=1.0)

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

