# ----------------------------------------------------------
# Library of cell classes and functions
#
# Tiago Branco, MRC Laboratory of Molecular Biology, 2013
# email: tbranco@mrc-lmb.cam.ac.uk
#
# Modified by 
# Balazs B Ujfalussy, MTA KOKI, 2018
# email: balazs.ujfalussy@gmail.com
#
# ----------------------------------------------------------

import numpy as np
import neuron

from neuron import h
from neuron import load_mechanisms
from neuron import gui

#load_mechanisms('/directory_where_mod_files_have_been_compiled')
h('objref nil')

# ----------------------------------------------------------
# MODELS
class L23(object):

    # Cell morphology is from cat, all lengths and diameters
    # are scaled to 70% to approximate it to mouse values

    def __init__(self):
        h('xopen("./L23.hoc")')
        props(self)
        self._geom()
        self._topol()
        self._changeLength()
        self._biophys()

    def _geom(self):
        self.axon = h.Section()
        # self.axon.L = 1
        self.axon.L = 300
        self.axon.diam = 1

    def _topol(self):            
        self.soma = h.soma
        self.dends = []
        for sec in h.allsec():
            self.dends.append(sec)
            n_seg = int(max(1, round(sec.L / self.seclength + 0.5)))
            if ((n_seg - int(n_seg/2)*2)==0) :
                n_seg = n_seg + 1
            sec.nseg = max(n_seg, 7) # minimum 7 segments to allow for hotspots
        self.dends.pop()   # Remove soma from the list
        self.dends.pop()   # and the Axon
        for sec in self.dends:
            sec.diam = sec.diam * self.rescale
        self.axon.connect(self.soma,1,0)
        
    
    def _biophys(self):
        for sec in h.allsec():
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA

    def _changeLength(self):
        for sec in h.allsec():
            sec.L = sec.L * self.rescale



############################################################################
## define the CA1 cell

class CA1(object):

    def __init__(self):
        h('xopen("./CA1.hoc")')
        propsCA1(self)
        # self._geom()
        self._topol()
        self._biophys()

    # def _geom(self):
    #     self.axon = h.Section()
    #     # self.axon.L = 1
    #     self.axon.L = 300
    #     self.axon.diam = 1

    def _topol(self):            
        self.soma = h.soma
        self.hill = h.hill
        self.iseg = h.iseg
        self.node = h.node
        self.inode = h.inode
        self.dends = []
        for sec in h.allsec():
            self.dends.append(sec)
            n_seg = int(max(1, round(sec.L / self.seclength + 0.5)))
            if ((n_seg - int(n_seg/2)*2)==0) :
                n_seg = n_seg + 1
            sec.nseg = n_seg
        for i in np.arange(8):
            self.dends.pop()   # Remove soma and axons from the list
        # self.axon.connect(self.soma,1,0)

    
    def _biophys(self):
        for sec in h.allsec():
            sec.cm = self.CM
            sec.insert('pas')
            sec.e_pas = self.E_PAS
            sec.g_pas = 1.0/self.RM
            sec.Ra = self.RA

        h.soma.g_pas = 1.0 /self.RM_soma 

        h.node[0].g_pas = 1.0 /self.RM_node 
        h.node[0].cm = self.CM

        h.node[1].g_pas = 1.0 /self.RM_node 
        h.node[1].cm = self.CM

        h.inode[0].g_pas = 1.0 /self.RM 
        h.inode[0].cm = self.CM_inode

        h.inode[1].g_pas = 1.0 /self.RM
        h.inode[1].cm = self.CM_inode

        h.inode[2].g_pas = 1.0 /self.RM
        h.inode[2].cm = self.CM_inode


        ## compensate for spines
        h('access soma')
        h('distance()')
        for sec in self.dends:
        # for sec in h.all_apicals: # in the Katz model compensation was done only for apicals.
            nseg = sec.nseg
            iseg = 0
            for seg in sec :
            # 0. calculate the distance from soma
                xx = iseg * 1.0/nseg + 1.0 / nseg / 2.0
                xdist=h.distance(xx, sec=sec)
                # print(sec.name(), xx, xdist)
                # 1. calculate the diameter of the segment
                xdiam=seg.diam
                # print(sec.name(), xx, xdiam)
                if ((xdist > self.spinelimit) + (xdiam < self.spinediamlimit)):
                    seg.cm = self.CM * self.spinefactor
                    seg.g_pas = self.spinefactor * 1.0/self.RM#_dend
                    # sec.Ra = self.RA_dend # Ra is a section variable...
                iseg = iseg + 1

############################################################################
## define the CA3 cell

def instantiate_swc(filename):
    """load an swc file and instantiate it"""
       
    # a helper library, included with NEURON
    h.load_file('import3d.hoc')
    
    # load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
    # Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
    # Import3d_Eutectic_read for Eutectic. (There is also an 
    # Import3d_Neurolucida_read for old Neurolucida files, but I've never seen one
    # in practice; try Import3d_Neurolucida3 first.)
    cell = h.Import3d_SWC_read()
    cell.input(filename)

    # easiest to instantiate by passing the loaded morphology to the Import3d_GUI
    # tool; with a second argument of 0, it won't display the GUI, but it will allow
    # use of the GUI's features
    i3d = h.Import3d_GUI(cell, 0)
    i3d.instantiate(None)

# ----------------------------------------------------------
# INSTRUMENTATION FUNCTIONS
def props(model):

    # morphology - rescale factor
    model.rescale = 0.7
    model.seclength = 10      # um, the length of a section
   
    # Passive properties
    model.CM = 1.0 
    model.RM = 7000.0 # 7000
    model.RA = 100.0 # 100
    model.E_PAS = -75
    model.CELSIUS = 35
    model.v_init = -75

    # Active properties
    model.Ek = -90
    model.Ena = 60
    model.Eca = 140
    
    model.gna_axon = 1000
    model.gkv_axon = 100
    
    model.gna_soma = 1000
    model.gkv_soma = 300 # 100
    model.gkm_soma = 2.2 
    model.gkca_soma = 3 
    model.gca_soma = 0.5 
    model.git_soma = 0.0003 
    
    model.gna_dend = 60 # 80
    model.gna_dend_hotSpot = 600 # 600
    model.gkv_dend = 3
    # model.gkv_dend_hotSpot = 60
    model.gkm_dend = 1
    model.gkca_dend = 3
    model.gca_dend = 0.5
    model.git_dend = 0.00015 
    model.gh_dend = 0


def propsCA1(model):
    
   # Passive properties
    model.CELSIUS = 35
    model.v_init = -68.3 # -68.3 for theta and -72 for replay
    model.RA = 100 # 150.00           # internal resistivity in ohm-cm 
    model.RA_dend = 100 # 200
    model.CM = 1 # 1                # 0.75     # specific membrane capacitance in uF/cm^2 
    model.CM_inode=0.04           # capacitance in myelin 

    model.RM = 20000 # was 3000 for hCond;  20000       # specific membrane resistivity at in ohm-cm^2 
    model.RM_dend = 20000 # 10 000 - 40 000; only used in Robustness analysis in Adam's paper
    model.RM_soma = 40000 # was 3000 for hCond;  20000       # specific membrane resistivity at the soma in ohm-cm^2 
    model.RM_inode = 40000 #200000          # inter-nodal resistivity with myelin
    model.RM_node = 50 #200000          # nodal resistivity

    model.E_PAS = -66 # -66 - set to v_init if passive       
    model.seclength = 10      # um, the length of a section
    model.spinefactor = 2       # 2 factor by which to change passive properties
    model.spinelimit = 100      # 100 distance beyond which to modify passive properties to account for spines 
    model.spinediamlimit = 1      # 100 distance beyond which to modify passive properties to account for spines 

   # Active properties - Values from the Spruston-lab (Katz et al., 2009) fitted only to trunk data!
    model.gna = 0.03      # 0.03; 0.01 sodium conductance in terminal branches
    model.gna_trunk = 0.04      # 0.04 sodium conductance
    model.gna_axon = 0.04      # 0.04 sodium conductance in the axon
    model.gna_soma = 0.2      # 0.04 - 0.2 sodium conductance in the soma
    model.gna_node = 15      # 30 - 15 sodium conductance in the axon
    model.nalimit = 500
    model.naslope = 0.001       # 0.001 is 'strong' propagation on the trunk
    model.gna_dend_hotSpot = 5

    model.gkdr = 0.02        # 0.005 delayed rectifier density in terminal branches 
    model.gkdr_trunk = 0.040        # 0.04 delayed rectifier density in the trunk
    model.gkdr_soma = 0.04        # 0.04 delayed rectifier density at the soma
    model.gkdr_axon = 0.04        # 0.04 delayed rectifier density at the axon

    model.gka = model.gkdr      # 0.005 A-type potassium density in terminal branches
    model.gka_trunk = 0.048      # 0.048  A-type potassium starting density in the trunk
    model.dlimit = 500        # cut-off for increase of A-type density 
    model.dprox = 100          # distance to switch from proximal to distal type 
    model.dslope=0.01         # slope of A-type density 


def init_active(model, axon=False, soma=False, dend=True, dendNa=False,
                dendCa=False):
    if axon:
        model.axon.insert('na'); model.axon.gbar_na = model.gna_axon
        model.axon.insert('kv'); model.axon.gbar_kv = model.gkv_axon
        model.axon.ena = model.Ena
        model.axon.ek = model.Ek
        print('active conductances added in the axon')
        
    if soma:
        model.soma.insert('na'); model.soma.gbar_na = model.gna_soma
        model.soma.insert('kv'); model.soma.gbar_kv = model.gkv_soma
        model.soma.insert('km'); model.soma.gbar_km = model.gkm_soma
        model.soma.insert('kca'); model.soma.gbar_kca = model.gkca_soma
        model.soma.insert('ca'); model.soma.gbar_ca = model.gca_soma
        model.soma.insert('it'); model.soma.gbar_it = model.git_soma
        model.soma.insert('cad');
        model.soma.ena = model.Ena
        model.soma.ek = model.Ek
        model.soma.eca = model.Eca
        print('somatic active conductances enabled')
        
    if dend:
        for d in model.dends:
            d.insert('na'); d.gbar_na = model.gna_dend*dendNa
            d.insert('kv'); d.gbar_kv = model.gkv_dend
            d.insert('km'); d.gbar_km = model.gkm_dend
            d.insert('kca'); d.gbar_kca = model.gkca_dend
            d.insert('ca'); d.gbar_ca = model.gca_dend*dendCa
            d.insert('it'); d.gbar_it = model.git_dend*dendCa
            d.insert('cad')
            d.ena = model.Ena
            d.ek = model.Ek
            d.eca = model.Eca


        print('active dendrites enabled', dendNa, model.gna_dend)

def init_passiveCA1(model):
    for sec in h.allsec():
        sec.e_pas = model.v_init

def init_activeCA1(model, soma=True, dend=True):

    if soma:
        model.soma.insert('nax'); model.soma.gbar_nax = model.gna_soma
        model.soma.insert('kdr'); model.soma.gkdrbar_kdr = model.gkdr_soma
        model.soma.insert('kap'); model.soma.gkabar_kap = model.gka

        model.hill.insert('nax'); model.hill.gbar_nax = model.gna_axon
        model.hill.insert('kdr'); model.hill.gkdrbar_kdr = model.gkdr_axon
        model.soma.insert('kap'); model.soma.gkabar_kap = model.gka

        model.iseg.insert('nax'); model.iseg.gbar_nax = model.gna_axon
        model.iseg.insert('kdr'); model.iseg.gkdrbar_kdr = model.gkdr_axon
        model.iseg.insert('kap'); model.soma.gkabar_kap = model.gka

        model.node[0].insert('nax'); model.node[0].gbar_nax = model.gna_node
        model.node[0].insert('kdr'); model.node[0].gkdrbar_kdr = model.gkdr_axon
        model.node[0].insert('kap'); model.node[0].gkabar_kap = model.gka*0.2

        model.node[1].insert('nax'); model.node[1].gbar_nax = model.gna_node
        model.node[1].insert('kdr'); model.node[1].gkdrbar_kdr = model.gkdr_axon
        model.node[1].insert('kap'); model.node[1].gkabar_kap = model.gka*0.2

        model.inode[0].insert('nax'); model.inode[0].gbar_nax = model.gna_axon
        model.inode[0].insert('kdr'); model.inode[0].gkdrbar_kdr = model.gkdr_axon
        model.inode[0].insert('kap'); model.inode[0].gkabar_kap = model.gka*0.2

        model.inode[1].insert('nax'); model.inode[1].gbar_nax = model.gna_axon
        model.inode[1].insert('kdr'); model.inode[1].gkdrbar_kdr = model.gkdr_axon
        model.inode[1].insert('kap'); model.inode[1].gkabar_kap = model.gka*0.2

        model.inode[2].insert('nax'); model.inode[2].gbar_nax = model.gna_axon
        model.inode[2].insert('kdr'); model.inode[2].gkdrbar_kdr = model.gkdr_axon
        model.inode[2].insert('kap'); model.inode[2].gkabar_kap = model.gka*0.2


        print('somatic and axonal active conductances enabled')
 
    if dend:

        for d in model.dends:
            d.insert('nad'); d.gbar_nad = model.gna
            d.insert('kdr'); d.gkdrbar_kdr = model.gkdr
            d.insert('kap'); d.gkabar_kap = 0
            d.insert('kad'); d.gkabar_kad = 0

        h('access soma')
        h('distance()')

        ## for the apicals: KA-type depends on distance
        ## density is as in terminal branches - independent of the distance
        for sec in h.all_apicals:
            nseg = sec.nseg
            iseg = 0
            for seg in sec:
                xx = iseg * 1.0/nseg + 1.0 / nseg / 2.0
                xdist=h.distance(xx, sec=sec)
                if (xdist > model.dprox):
                    seg.gkabar_kad = model.gka
                else:
                    seg.gkabar_kap = model.gka
                iseg = iseg + 1
       
        h('access soma')
        h('distance()')

        ## distance dependent A-channel densities in apical trunk dendrites
        ##      1. densities increase till 'dlimit' with dslope
        ##      2. proximal channels switch to distal at 'dprox'
        ##      3. sodiom channel density also increases with distance
        for sec in h.primary_apical_list:
            nseg = sec.nseg
            sec.insert('nax')
            iseg = 0
            for seg in sec :
            # 0. calculate the distance from soma
                xx = iseg * 1.0/nseg + 1.0 / nseg / 2.0
                xdist=h.distance(xx, sec=sec)
            # 1. densities increase till 'dlimit' with dslope
                if (xdist > model.dlimit):
                    xdist = model.dlimit
            # 2. proximal channels switch to distal at 'dprox'
                if (xdist > model.dprox):
                    seg.gkabar_kad = model.gka_trunk*(1+xdist*model.dslope)
                else:
                    seg.gkabar_kap = model.gka_trunk*(1+xdist*model.dslope)
                iseg = iseg + 1
            # 3. sodiom channel density also increases with distance
                if (xdist>model.nalimit):
                    xdist=model.nalimit
                seg.gbar_nax = model.gna_trunk*(1+xdist*model.naslope)
                # print(sec.name(), model.gna_trunk, xdist, model.naslope, seg.gbar_nax)
                seg.gbar_nad = 0
                seg.gkdrbar_kdr = model.gkdr_trunk
 
        ## for the basals: all express proximal KA-type
        ## density does not increase with the distance
        for sec in h.all_basals:
            for seg in sec:
                seg.gkabar_kap = model.gka


        # for sec in h.parent_list:
        #     for seg in sec:
        #         seg.gbar_nad = model.gna * 0 
        #         sec.insert('nax')
        #         seg.gbar_nax = 30 * model.gna_trunk*(1+model.nalimit*model.naslope)

        print('active dendrites enabled')

def CA1_hotSpot(model, iden, gNa_factor=10):
    n_seg = model.dends[iden].nseg
    model.dends[iden].nseg = max(n_seg, 10)
    nmax = model.dends[iden].nseg

    # d = model.dends[iden]
    # d.insert('nafast'); d.gbar_nafast = 0

    s = 0
    min_spot = np.floor(2.*nmax/6.)
    max_spot = np.floor(4.*nmax/6.)
    for seg in model.dends[iden]:
        s+=1
        if ((s >= min_spot) & (s < max_spot)):
            seg.gbar_nad = model.gna * gNa_factor 
            # seg.gbar_nafast = model.gna * gNa_factor 
            # print s
    print('hotspots added')


def hotSpot(model):
    for section in model.dends:
        s = 0
        nmax = section.nseg
        min_spot = np.ceil(3.*nmax/7.)
        max_spot = np.ceil(4.*nmax/7.)
        for seg in section:
            if ((s > min_spot) & (s <= max_spot)):
                seg.gbar_na = model.gna_dend_hotSpot    
                # seg.gbar_kv = model.gkv_dend_hotSpot    
                # print section.name(), 'hotspot added', s
            else:
                seg.gbar_na = 0
                # seg.gbar_kv = 0
            s+=1
    print('hotspots added')


def add_somaStim(model, p=0.5, onset=20, dur=1, amp=0):
    model.stim = h.IClamp(model.soma(p))
    model.stim.delay = onset
    model.stim.dur = dur
    model.stim.amp = amp    # nA

def add_dendStim(model, p=0.5, dend=10, onset=20, dur=1, amp=0):
    model.stim = h.IClamp(model.dends[dend](p))
    model.stim.delay = onset
    model.stim.dur = dur
    model.stim.amp = amp    # nA

def add_dendStim4(model, dends, onset=20, dur=1, amp=0):
    model.stim1 = h.IClamp(model.dends[dends[0]](0.5))
    model.stim1.delay = onset
    model.stim1.dur = dur
    model.stim1.amp = (1) * amp   # nA

    model.stim2 = h.IClamp(model.dends[dends[1]](0.5))
    model.stim2.delay = onset
    model.stim2.dur = dur
    model.stim2.amp = (1) * amp   # nA

    model.stim3 = h.IClamp(model.dends[dends[2]](0.5))
    model.stim3.delay = onset
    model.stim3.dur = dur
    model.stim3.amp = (1) * amp   # nA

    model.stim4 = h.IClamp(model.dends[dends[3]](0.5))
    model.stim4.delay = onset
    model.stim4.dur = dur
    model.stim4.amp = (1) * amp    # nA

def synDist(model,locs):
    nsyn = len(locs)
    DSyn = np.zeros([nsyn, nsyn])
    fromSyn = 0
    for loc in locs:
        fromDend = loc[0]
        fromX = loc[1]
        fromSection = model.dends[fromDend]
        h.distance(0, fromX, sec=fromSection)

        toSyn = 0
        for toLoc in locs:
            toDend = toLoc[0]
            toX = toLoc[1]
            toSection = model.dends[toDend]
            x = h.distance(toX, sec=toSection)
            DSyn[toSyn, fromSyn] = x
            toSyn = toSyn + 1
        fromSyn = fromSyn + 1
    return DSyn

def surface_area(model):
    sa = model.soma.diam * model.soma.L * np.pi
    ndends = len(model.dends)
    for i in range(ndends):
        sa = sa + model.dends[i].diam * model.dends[i].L * np.pi
        print('dendrite', model.dends[i].name(), 'diam', model.dends[i].diam, 'length', model.dends[i].L)
    return sa



def add_syns(model, data):
    print('adding synapses using the new function!')
    model.AMPAlist = []
    model.ncAMPAlist = []
    AMPA_gmax = data.Agmax/1000.   # Set in nS and convert to muS

    if (data.SPINES):
        # h.nspines = len(data.Elocs)
        h('nspines = 0')
        h.nspines = len(data.Elocs)
        print('we will create', h.nspines, 'dendritic spines for excitatory synapses')
        h('create shead[nspines]')
        h('create sneck[nspines]')

    if (data.NMDA):
        model.NMDAlist = []
        model.ncNMDAlist = []
        NMDA_gmax = data.Ngmax/1000.   # Set in nS and convert to muS

    spi=0 # indexing the spines - in hoc

    for loc in data.Elocs:
        locInd = int(loc[0])
        if (locInd == -1):
            synloc = model.soma
        else:
            if (data.SPINES):

                neck = h.sneck[spi]
                neck.L = data.sneck_len
                neck.diam = data.sneck_diam
                neck.insert("pas")
                neck.e_pas = model.E_PAS
                neck.g_pas = 1.0/model.RM
                neck.Ra = model.RA
                neck.cm = model.CM

                head = h.shead[spi]
                head.L = data.shead_len
                head.diam = data.shead_diam
                head.insert("pas")
                head.e_pas = model.E_PAS
                head.g_pas = 1.0/model.RM
                head.Ra = model.RA
                head.cm = model.CM

                head.connect(neck, 1, 0)
                neck.connect(model.dends[int(loc[0])], loc[1], 0)
                synloc = h.shead[spi]
                synpos = 0.5
                spi = spi + 1
            else:
                synloc = model.dends[int(loc[0])]
                synpos = float(loc[1])
        # print loc[0], loc[1]
        AMPA = h.Exp2Syn(synpos, sec=synloc)
        AMPA.tau1 = data.Atau1
        AMPA.tau2 = data.Atau2
        NC = h.NetCon(h.nil, AMPA, 0, 0, AMPA_gmax) # NetCon(source, target, threshold, delay, weight)
        model.AMPAlist.append(AMPA)
        model.ncAMPAlist.append(NC)

        if (data.NMDA):
            NMDA = h.Exp2SynNMDA(synpos, sec=synloc) 
            NMDA.tau1 = data.Ntau1
            NMDA.tau2 = data.Ntau2
            NC = h.NetCon(h.nil, NMDA, 0, 0, NMDA_gmax)
            x = float(loc[1])
            model.NMDAlist.append(NMDA)
            model.ncNMDAlist.append(NC)   


    print('AMPA synapses added')
    if (data.NMDA):
        print('dExp NMDA synapses generated')

    if (data.GABA):
        model.GABAlist = []
        model.ncGABAlist = []
        GABA_gmax = data.Igmax/1000.   # Set in nS and convert to muS
    
        if (data.GABA_B):
            model.GABA_Blist = []
            model.ncGABA_Blist = []
            GABAB_gmax = data.Bgmax/1000.   # Set in nS and convert to muS

        for loc in data.Ilocs:
            locInd = int(loc[0])
            if (locInd == -1):
                synloc = model.soma
            else:
                synloc = model.dends[int(loc[0])]
            GABA = h.Exp2Syn(float(loc[1]), sec=synloc) 
            GABA.tau1 = data.Itau1
            GABA.tau2 = data.Itau2
            GABA.e = data.Irev
            NC = h.NetCon(h.nil, GABA, 0, 0, GABA_gmax)
            model.GABAlist.append(GABA)
            model.ncGABAlist.append(NC)

            if (data.GABA_B):
                GABAB = h.Exp2Syn(float(loc[1]), sec=synloc) 
                GABAB.tau1 = data.Btau1
                GABAB.tau2 = data.Btau2
                GABAB.e = data.Brev
                NC = h.NetCon(h.nil, GABAB, 0, 0, GABAB_gmax)
                model.GABA_Blist.append(GABAB)
                model.ncGABA_Blist.append(NC)
    
        print('inhibitory synapses generated')
        if (data.GABA_B):
            print('GABA_B synapses generated')

# ----------------------------------------------------------
# SIMULATION RUN
def simulate(model, t_stop=100, recDend=False, i_recDend=11, x_recDend=0.5):
    trec, vrec = h.Vector(), h.Vector()
    vDendRec = []
    trec.record(h._ref_t)
    vrec.record(model.soma(0.5)._ref_v)

    if recDend:
        n = 0
        for i_dend in i_recDend:
            vDendRec.append(h.Vector())
            x_dend = x_recDend[n]
            vDendRec[n].record(model.dends[i_dend](x_dend)._ref_v)
            n+=1
         

    h.celsius = model.CELSIUS
    h.finitialize(model.v_init)
    neuron.run(t_stop)
    return np.array(trec), np.array(vrec), np.array(vDendRec)
