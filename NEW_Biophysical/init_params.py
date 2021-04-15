# Simulation control parameters
data.dt = 0.2
lb.h.dt = data.dt
lb.h.steps_per_ms = 1.0/lb.h.dt

#----------------------------------------------------------------------------
### main experemient parameters
data.model = 'CA1' # L23 or CA1
data.stimType = 'place' # poisson, place, replay, minis, nIter, SStim, DStim
data.synType = 'clust2' # synapse distribution - NULL, single, alltree, clustered, clust2, local_regular, global_regular_lrand, global_regular_lreg
data.recordDend = True

data.placeType = 'random_NR' # balanced, random_N (number of neurons), random_NR (number and peak firing rate)
data.randomW = False # True or False - weights are coming from a random distribution
if (data.placeType == 'balanced'): data.randomW = False


#----------------------------------------------------------------------------
# cell type parameters: active or passive?

data.ACTIVE = False
data.ACTIVE_soma = False
data.ACTIVE_dend = False
data.ACTIVEhotSpot = False
if (data.actType == 'aDend'):
    data.ACTIVE = True
    data.ACTIVE_dend = True
    data.ACTIVEhotSpot = data.AHS
if (data.actType == 'aSoma'):
    data.ACTIVE = True
    data.ACTIVE_soma = True
if (data.actType == 'active'):
    data.ACTIVE = True
    data.ACTIVE_dend = True
    data.ACTIVE_soma = True
    data.ACTIVEhotSpot = data.AHS

data.hotspot_branches = [5, 6, 52, 53, 56, 57, 59, 60, 82, 87, 111, 113]

#----------------------------------------------------------------------------
## FIRST: set stimulation parameters for minis or Niter
#----------------------------------------------------------------------------
data.st_onset = 0.1 # in seconds - only for minis
nseg_den = 10 # except when nIter

if (data.stimType == 'nIter'):
    if (data.synType != 'single'):
        print('error: synapse distribution must be single for nIter!')
        sys.exit(1)

    data.tInterval = .3
    data.TSTOP = 0.3
    data.GABA = False
    nseg_den = 50
    data.nIter = 1 # 2 



if (data.stimType == 'minis'):
    data.TSTOP = 0.35
    data.GABA = True


#----------------------------------------------------------------------------
# SECOND: input parameters for place cells stimulus
#----------------------------------------------------------------------------
data.NMDA = True
data.GABA = True   
data.GABA_B = True

data.SPINES = False
data.Lmin = 60
weight_factor_A = 1 # multiplcative weight scale of the clustered synapses - AMPA
weight_factor_N = 1 # multiplcative weight scale of the clustered synapses - NMDA
data.g_factor = 1 # all synapses are scaled with this factor


data.Ensyn = 2000
data.Insyn = 200
data.Erate = 0.5
data.Irate = 7.4

if (data.stimType == 'minis'):
    data.Ensyn = 100
    data.Insyn = 10

if (data.synType == 'single'):
    # synapse location
    data.Ensyn = 30
    ## this part of the branch will be covered by synapses in a single branch stimulation experiment
    data.locSeg = [0, 1]


data.cdends = [[8, 13, 36, 59], [5, 56, 144, 147], [20, 21, 45, 109], [38, 111, 127, 137], [6, 96, 146, 147], [5, 36, 120, 137], [28, 56, 57, 137], [9, 48, 82, 92], [28, 96, 98, 100], [32, 93, 108, 109]]
# data.cdends = [[5, 56, 82, 111], [5, 56, 144, 147], [20, 21, 45, 109], [38, 111, 127, 137], [6, 96, 146, 147], [5, 36, 120, 137], [28, 56, 57, 137], [9, 48, 82, 92], [28, 96, 98, 100], [32, 93, 108, 109]]
# if data.ACTIVEhotSpot: 
#     data.cdends = [data.hotspot_branches]

data.locDend = data.cdends[PAR2] # location of the clusters


# trunk dendrites: [-1, 76, 80, 88, 90, 94, 102, 104, 106, 110]
# their lengths:   [57, 57, 21, 7,  31, 11, 33,  7,   5,   20]
data.locDendRec = data.locDend
data.xDendRec = [0.5] * (len(data.locDend) + 1)


#----------------------------------------------------------------------------
# THIRD: synapse parameters
#----------------------------------------------------------------------------

data.Agmax = 0.6 * data.g_factor# 0.1 nS
data.Atau1 = 0.1 # ms
data.Atau2 = 1 # ms

data.Ngmax = 0.8 * data.g_factor # 0.8 nS
data.Ntau1 = 2 # ms
data.Ntau2 = 50 # ms

data.Itau1 = 0.1   # ms
data.Itau2 = 4  # ms
data.Irev = -65 # -65
data.Igmax = 0.1 * data.g_factor   # 0.7 nS

data.Btau1 = 1   # ms
data.Btau2 = 40   # ms
data.Brev = -80 # -75
data.Bgmax = 0.1 * data.g_factor  # 1.2 nS

if (data.g_factor == 2):
    data.Igmax = 0.85 * data.g_factor   # 0.85 nS
    data.Brev = -75 # -75
    data.Bgmax = 1.45 * data.g_factor  # 1.45 nS

#----------------------------------------------------------------------------
# FOURTH: spines
#----------------------------------------------------------------------------
# Spine dimesions - spines are hoc objects, so these parameters will be needed in the hoc environment
# data.sneck_diam = 0.04 # um (high resistance), approximate Grunditz et al model */
# sneck_diam = 0.15 # um (low resistance) */
data.sneck_diam = 0.077 # um  from Harnett et al, 2012
data.sneck_len = 1.58 # 1.58 um
data.shead_diam = 0.5 # um
data.shead_len = 0.5 # um

## decreasing the amount of inhibition by removing inhibitory spikes
data.removeIspikes = 0 # 0-1. The probability of removing an inhibitory spike