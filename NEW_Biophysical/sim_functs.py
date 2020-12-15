#----------------------------------------------------------------------------
# Functions and Classes
import struct
import sys
from numpy import loadtxt
import copy

def initOnsetSpikes():
    model.ncAMPAlist[0].event(data.st_onset * 1000)

def initSpikes():
    if (len(data.etimes)>0):
        for s in data.etimes:
            model.ncAMPAlist[int(s[0])].event(float(s[1]))
            if data.NMDA: model.ncNMDAlist[int(s[0])].event(float(s[1]))
    if data.GABA == True:
        if (len(data.itimes)>0):
            # print data.itimes[1]
            for s in data.itimes:
                model.ncGABAlist[int(s[0])].event(float(s[1]))
                if data.GABA_B: model.ncGABA_Blist[int(s[0])].event(float(s[1]))

def storeSimOutput(v,vD):
        data.vdata.append(v)
        data.vDdata.append(vD)

def storeSimInputOutput(v,vD,Et,It):
        data.vdata.append(v)
        data.vDdata.append(vD)
        if (len(data.stim) == 0):
            nrep = 1
        else :
            nrep = max(data.stim[:,0]) + 1
            
        if not Et is None:
            Et = np.column_stack((nrep * np.ones(len(Et[:,0])), Et))
            if (len(data.stim) == 0):
                data.stim = Et
            else :
                data.stim = np.row_stack((data.stim, Et))
        if not It is None :
            It = np.column_stack((-1 * nrep * np.ones(len(It[:,0])), It))
            if (len(data.stim) == 0):
                data.stim = It
            else :
                data.stim = np.row_stack((data.stim, It))


#-----------------------------------------------
# Synapse location functions
#-----------------------------------------------

def genDendLocs(dends, nsyn, spread):
    # insert nsyn synapses to dendrites dends, uniform spread within a branch
    locs = []
    n_dends = len(dends)
    if isinstance(nsyn, list):
        nsyn = np.repeat(nsyn[0], len(dends))
    else :
        nsyn = np.repeat(nsyn, len(dends))        
    for i_dend in np.arange(0,n_dends):
        dend = dends[i_dend]
        nsyn_dend = nsyn[i_dend]
        isd = (spread[1]-spread[0])/float(nsyn_dend)
        pos = np.arange(spread[0], spread[1], isd)[0:nsyn_dend] 

        if (len(pos) != nsyn_dend):
            print('error: synapse number mismatch, stop simulation! dend:', i_dend, 'created=', len(pos), '!=', nsyn_dend)
            sys.exit(1)
        for p in pos:
            locs.append([dend, p])

    return locs

def genRandomLocs(nsyn, seed, dend_ids=None):
    # randomly choose a dendritic branch, proportionally to its length
    np.random.seed(seed)
    nden = len(model.dends)
    Ldends = np.zeros(nden)
    for ii in np.arange(0, nden): Ldends[ii] = model.dends[ii].L

    # bfname='Dend_length.bin'
    # binfile = file(bfname, 'wb')
    # # and write out two integers with the row and column dimension
    # header = struct.pack('2I', Ldends.shape[0], 1)
    # binfile.write(header)
    # ddata = struct.pack('%id' % Ldends.shape[0], *Ldends)
    # binfile.write(ddata)
    # binfile.close()

    # print (dend_ids)
    Pdends = Ldends / sum(Ldends)
    if dend_ids is not None:
        Pdends = np.zeros(nden)
        for dend in dend_ids:
            Pdends[dend] = Ldends[dend] / sum(Ldends)
        Pdends = Pdends / sum(Pdends)

    # print(Pdends)    
    Ndends = np.random.multinomial(nsyn, Pdends)

    #  and insert a nsyn synapses at random locations within the branch
    locs = []
    for dend in np.arange(0, nden): # 
        nsyn_dend = Ndends[dend]
        if (nsyn_dend > 0):
            locs_dend = np.sort(np.random.uniform(0, 1, nsyn_dend))
            for s in np.arange(0,nsyn_dend):
                locs.append([dend, locs_dend[s]])
    # print locs
    return locs

def genClusts(Nclust, Ncell_per_clust, minL, seed, clocs=None):
    # randomly choose a dendritic branch, larger than minL for clustered synapses
    # clocs: 
    if (minL < Ncell_per_clust):
        print('error: cluster size mismatch, stop simulation! Number of cells per cluster (1 um inter-spine distance):', Ncell_per_clust, ', minL:', minL)
        sys.exit(1)
    np.random.seed(seed)
    nden = len(model.dends)
    Ldends = np.zeros(nden)
    for ii in np.arange(0, nden): Ldends[ii] = model.dends[ii].L
    # idends = np.flatnonzero((Ldends < minL) & (Ldends > 10)) # to put clusters on short segments
    idends = np.flatnonzero(Ldends > minL)
    # print(idends)

    replace_clusts = False
    if (Nclust > len(idends)):
        print('warning: cluster number mismatch, multiple clusters are allocated to the same branch! Nclust:', Nclust, ', number of available branches:', len(idends), ', minL:', minL)
        replace_clusts = True
    
    # if some dendrites are preselected ...
    if clocs is not None:
        k = len(clocs)
        clDends1 = clocs # we choose those dendrites
        if (Nclust > k) : # we choose randomly from the others
            if (replace_clusts == False):
                idends = np.setdiff1d(idends, clocs)
            clDends2 = np.random.choice(idends, Nclust-k, replace=replace_clusts)
            clDends = np.concatenate((clDends1, clDends2))
        else:
            clDends = np.array(clDends1)

    else:
        clDends = np.random.choice(idends, Nclust, replace=replace_clusts)
    #  and insert a nsyn synapses at random locations within the branch
    # print(clDends)
    locs = []
    i_dend = 0
    for dend in clDends: # 
        locstart = np.random.uniform(0, 1 - Ncell_per_clust / Ldends[dend], 1)
        if clocs is not None:
            if i_dend < k:
                locstart = (1 - Ncell_per_clust / Ldends[dend])/2
        for s in np.arange(0,Ncell_per_clust):
            locs.append([dend, np.asscalar(s/Ldends[dend] + locstart)]) # 1 um distance between spines
        i_dend = i_dend + 1
    return locs

##########################
## these functions are needed to distribute cells assemblies uniformly within a branch
## for L synapses on a branch first select a number P that is
##      - co-prime with L
##      - somehwat smaller than L/4
## the idea is to construct a {sequence with increments of P} %% L
## and use it to index the contacts
## we only do this if L > 6

def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

def genClustLocs(nsyn, Nclust, Ncell_per_clust, seed, midle=False):
    # Nclust clusters and random background innervation
    # cluster locations are SELECTED from a random background
    # clustering: 
    #     1. select pre clusters by genClustStarts
    #     2. select post cluster locations by genClustStarts
    #     3. generate ordered Elocs
    #     4. remove post cluster locations from Elocs
    #     5. cycle through all clusters and add the post from either the random or the post cluster locations
    # output: Elocs: list with [dend id, location]

    # 1. presynaptic partners belonging to the clusters
    if midle == True :
        Nsyn_in_clust = Nclust * Ncell_per_clust# number of synapses in clusters
        Nsyn_rand = nsyn - Nsyn_in_clust
        istart = int(Nsyn_rand / 2)
        iend = int(Nsyn_rand / 2) + Nsyn_in_clust
        prestarts = np.arange(istart, iend, Ncell_per_clust)
    else:
        prestarts = genClustStarts(nsyn, Nclust, Ncell_per_clust, 100 * seed)

    # index of pre cells IN clusters
    ind_clustpre = np.arange(prestarts[0], prestarts[0]+Ncell_per_clust)
    if Nclust > 0:
        for j in np.arange(1, Nclust):
            istart = prestarts[j]
            iend = prestarts[j] + Ncell_per_clust
            ind_clustpre = np.concatenate((ind_clustpre, np.arange(istart, iend)))
    
    # print ind_clustpre

    # 2. synapses belonging to the clusters
    poststarts = genClustStarts(nsyn, Nclust, Ncell_per_clust, 100 * seed + 1)

    # index of synapses NOT in clusters
    ind_randpost = np.arange(poststarts[0])
    if Nclust > 1:
        for j in np.arange(1, Nclust):
            istart = poststarts[j-1] + Ncell_per_clust
            iend = poststarts[j]
            ind_randpost = np.concatenate((ind_randpost, np.arange(istart, iend)))
    else :
        j = 0
    istart = poststarts[j] + Ncell_per_clust
    iend = nsyn
    ind_randpost = np.concatenate((ind_randpost, np.arange(istart, iend)))
    print(len(ind_randpost))

    # shuffling the post clusters
    np.random.seed(100*seed + 2)
    np.random.shuffle(poststarts)

    # index of synapses IN clusters
    ind_clustpost = np.arange(poststarts[0], poststarts[0]+Ncell_per_clust)
    if Nclust > 1:
        for j in np.arange(1, Nclust):
            istart = poststarts[j]
            iend = poststarts[j] + Ncell_per_clust
            ind_clustpost = np.concatenate((ind_clustpost, np.arange(istart, iend)))
    print(len(ind_clustpost))

    # 3. random locations - all synapses
    ordered_Elocs = genRandomLocs(nsyn, seed=100*seed+3)

    # 4. random locations - NOT in clusters
    random_Elocs = [ ordered_Elocs[i] for i in ind_randpost]
    np.random.seed(100*seed + 4)
    np.random.shuffle(random_Elocs)

    # random locations - IN clusters
    clustered_Elocs = [ ordered_Elocs[i] for i in ind_clustpost]
    
    # 5. cycle through all clusters and add the post from either the 
    #       random or the post cluster locations
    j = 0 # POST - clustered
    k = 0 # POST - random
    Elocs = ordered_Elocs
    for i in np.arange(nsyn): # PRE
        if i in ind_clustpre: # take it from the cluster - list
            Elocs[i] = clustered_Elocs[j]
            j = j + 1
        else :
            Elocs[i] = random_Elocs[k]
            k = k + 1
    print(len(Elocs))

    return Elocs, ind_clustpre

def addClustLocs(nsyn, Nclust, Ncell_per_clust, seed, midle=False, clocs=None, Lmin=60):
    # Nclust clusters in a random background innervation
    # clusters are ADDED to the random background
    # clustering: 
    #     1. select pre clusters by genClustStarts
    #     2. generate random Elocs
    #     3. select post cluster locations by genClusts
    #     4. add post cluster locations
    # output: Elocs: list with [dend id, location]
    
    # 1. presynaptic partners belonging to the clusters
    if midle == True : # clusters start at the middle of the maze
        Nsyn_in_clust = Nclust * Ncell_per_clust# number of synapses in clusters
        Nsyn_rand = nsyn - Nsyn_in_clust
        istart = 880 #int(Nsyn_rand / 2)
        iend = 1120 # int(Nsyn_rand / 2) + Nsyn_in_clust
        # istart = int(Nsyn_rand / 2)
        # iend = int(Nsyn_rand / 2) + Nsyn_in_clust
        prestarts = np.arange(istart, iend, Ncell_per_clust)        
    else:
        prestarts = genClustStarts(nsyn, Nclust, Ncell_per_clust, 100 * seed)

    # print(prestarts)
    # index of pre cells IN clusters
    ind_clustpre = np.arange(prestarts[0], prestarts[0]+Ncell_per_clust)
    if Nclust > 0:
        for j in np.arange(1, Nclust):
            istart = prestarts[j]
            iend = prestarts[j] + Ncell_per_clust
            ind_clustpre = np.concatenate((ind_clustpre, np.arange(istart, iend)))

    # 2. random locations - background synapses
    if (Nsyn_rand): 
        bg_Elocs = genRandomLocs(Nsyn_rand, seed=100*seed+1)
        np.random.seed(100*seed + 2)
        np.random.shuffle(bg_Elocs)
    else:
        bg_Elocs = []

    # 3. clustered locations
    # def genClusts(Nclust, Ncell_per_clust, minL, seed):
    clustered_Elocs = genClusts(Nclust, Ncell_per_clust, Lmin, seed=100*seed+4, clocs=clocs)
    # clustered_Elocs = genClusts(Nclust, Ncell_per_clust, 1, seed=100*seed+4)
    
    # 4. cycle through all clusters and add the post from either the 
    #       random or the post cluster locations
    j = 0 # POST - clustered
    k = 0 # POST - random
    Elocs = bg_Elocs + clustered_Elocs
    for i in np.arange(nsyn): # PRE
        if i in ind_clustpre: # take it from the cluster - list
            Elocs[i] = clustered_Elocs[j]
            j = j + 1
        else :
            Elocs[i] = bg_Elocs[k]
            k = k + 1
    print(len(Elocs))

    return Elocs, ind_clustpre




#-----------------------------------------------
# Input generation functions
#-----------------------------------------------

def readTrain(type, Ensyn, Insyn, Erate, Irate, duration, stimseed, rep):

    ## eitimes = readPoissonTrain(duration, ir, rep, Ensyn, bgErate, fgErate, Insyn, bgIrate, fgIrate)
    ## rates are in Hz
    ## duration: total simulation time in s
    ## times: two columns: cell id, spike times (ms)
    Etimes = np.array([])
    print('read Train: duration: ', duration)

    fname = './'+str(data.model)+'stims/'+str(type)+'/Espikes_d'+str(duration)+'_Ne'+str(Ensyn)+'_Re'+str(Erate)+'_rseed'+str(stimseed)+'_rep'+str(rep)+'.dat'
    print('loading input from file', fname)
    Etimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
    
    if (Insyn > 0) :
        Itimes = np.array([])
        fname = './'+str(data.model)+'stims/'+str(type)+'/Ispikes_d'+str(duration)+'_Ni'+str(Insyn)+'_Ri'+str(Irate)+'_rseed'+str(stimseed)+'_rep'+str(rep)+'.dat'
        Itimes = loadtxt(fname, comments="#", delimiter=" ", unpack=False)    
        times = [Etimes, Itimes]
    else :
        times = Etimes

    return times



#--------------------------------------------------
# Simulation functions
#--------------------------------------------------


def sim_PlaceInput(type, Ensyn, Insyn, Erate, Irate, duration, stimseed, rep, elimIspike=0):
    ## spike train read from file - Erate [Irate] is the rate of the individual input neuronss

    print('sim Place Input: duration: ', duration)
    eitimes = readTrain(type, Ensyn, Insyn, Erate, Irate, duration, stimseed, rep)
    if ((Insyn > 0) & (Irate > 0)) :
        etimes = eitimes[0]
        itimes = eitimes[1]
        #### this is a hack to phase shift inhibitory spikes
        # itimes[:,1] = (itimes[:,1] + 115) % max(itimes[:,1])
        print(len(etimes[:,0]), 'E and ',  len(itimes[:,0]),  'I spikes read from file')
        if (elimIspike > 0):
            N_ispikes = len(itimes)
            i_index = np.sort(np.random.choice(N_ispikes, int(round((1-elimIspike) * N_ispikes)), replace=False))
            itimes = itimes[i_index]
            print('number of I spikes reduced to', len(itimes))
        data.itimes = itimes
    else :
        etimes = eitimes
        print(len(etimes[:,0]), 'E spikes generated')
        itimes = None
        data.itimes = itimes

    data.etimes = etimes
 
    # Run
    fih = lb.h.FInitializeHandler(1, initSpikes)
    taxis, v, vD= lb.simulate(model, t_stop=data.TSTOP * 1000, recDend=data.recordDend, i_recDend=data.locDendRec, x_recDend=data.xDendRec)
    
    return taxis, v, vD, etimes, itimes


#--------------------------------------------------
# Iteration functions
#--------------------------------------------------


def SIM_PlaceIteration(type, Ensyn, Insyn, Erate, Irate, duration, Nrep, stimseed, elimIspike=1):
    # Nrep simulation with same paramaters but different input spikes read from file
    # rates are in Hz, duration in s
    
    for iter in np.arange(Nrep):
    # iter = 11
        print('Running iteration', iter, 'duration: ', duration)
        data.taxis, v, vD, Et, It = sim_PlaceInput(type=type, Ensyn=Ensyn, Insyn=Insyn, Erate=Erate, Irate=Irate, duration=duration, stimseed=stimseed, rep=iter, elimIspike=elimIspike)
        storeSimInputOutput(v,vD,Et,It)



def SIM_minisDistribution(data):
#   ---------------------------------------
#   first each excitatory synapse is activated alone
    It, Et, data.etimes, data.itimes = [], [], [], []
    Enum = data.Ensyn
    for nsyn in np.arange(0, Enum):
    # for nsyn in [23, 952, 981]:
        print('activating excitatory synapse ', nsyn, ' alone')
        data.etimes = np.array([[nsyn, data.st_onset * 1000]])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000, recDend=data.recordDend, i_recDend=data.locDendRec, x_recDend=data.xDendRec)
        Et = data.etimes
        storeSimOutput(v, vD)

#   ---------------------------------------
#   next, each inhibitory synapse is activated alone
    It, Et, data.etimes, data.itimes = [], [], [], []
    Inum = data.Insyn
    for nsyn in np.arange(0, Inum):
        print('activating inhibitory synapse ', nsyn, ' alone')
        data.itimes = np.array([[nsyn, data.st_onset * 1000]])
        fih = lb.h.FInitializeHandler(1, initSpikes)
        taxis, v, vD = lb.simulate(model, t_stop=data.TSTOP * 1000, recDend=data.recordDend, i_recDend=data.locDendRec, x_recDend=data.xDendRec)
        It = data.itimes
        storeSimOutput(v, vD)

    data.taxis = taxis
