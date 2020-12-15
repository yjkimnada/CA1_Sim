### General information

Library for simulating the effect of synaptic clustering on neuronal selectivity

for more details, see our paper on bioRxive:
**Balazs B Ujfalussy and Judit K Makara:** *Impact of functional synapse clusters on neuronal response selectivity* `https://doi.org/ 10.1101/634220`

author: Balazs B Ujfalussy, 2019 `balazs.ujfalussy@gmail.com`

This code reproduces the main findings of the paper: the sensitivity of CA1 and L23 neurons to the cluster size under 3 different input conditions (theta, sharp wave for CA1 and oriented grating stimulus for L23). The library contains only a limited set of in vivo-like stimuli (usually 1 or two random seeds). 

The scripts used for generating the stimuli can be found here:
https://bitbucket.org/bbu20/popact/

### Requirements

To run the code, you need python 2.7 installed with the usual sciantific packages as well as neuron installed as a python module. Follow the instructions on Neuron's website: `https://neuron.yale.edu/neuron/download/compilestd_osx`. Installation takes a few minutes if you already have python and XCode installed and updated. 

The code has been tested with Python 2.7.16, Neuron 7.4. installed on a Macbook Pro Mid 2015, macOS Mojave 10.14.5, 

To run the code on OsX you first compile the mod files with nrnivmodl ./ in the directory, then start python (ipython) and type `import mainAll`

### Details

in the file mainAll.py you can control the major properties of the simulation:
 
  * data.model: cell type, `L23` or `CA1`
 
  * data.stimType: stimulus type. Possible values are: 
	* `SStim`: somatic current clamp
    * `DStim`: dendritic current clamp
    * `nIter`: stimulating an increasing number of synapses on a single branch
    * `minis`: one synapse at a time to record the mini's distribution
    * `place`: place field or orientation selective input - read from files in the CA1stims or L23stims directory
    * `replay`: SPW inputs - read from files in the CA1stims directory; only if model == CA1
    * `poisson`: random Poisson spike train
 
  * data.synType: synapse distribution. Possible values:
    * `NULL` - for current clamp
    * `single`: single branch stimulation
    * `alltree`: random synapses distributed along the whole dendritic tree
    * `clustered`: clustered synapses + random background. Inter-synapse distance is similar between synapses within clusters and background synapses
    * `clust2`: ISD (inter/synapse distance) is set separately for clustered synapses
    * `clust3`: used only with SPWs, input place cells and background cells are separated
    * `local_regular`: Locally balanced synapses
    * `global_regular_lrand`: Globally balanced synapses, randomly reaaranged within branches
    * `global_regular_lreg`: Fully balanced synapses
 
  * data.actType: whether the neuron is active or passive. possible values:
    * `passive`: fully passive
    * `aSoma`: only the soma is active, dendrites are passive
    * `aDend`: only the dendrites are active, soma is passive
    * `active`: fully active
 
 Further parameters can be set in the file init_params.py. Typical run time is a few (1-2) minutes.
 
 Analysis of the produced data is implemented by the R scripts provided in the analysis folder (see the Readme file there)

### Expected output

 The results are saved into a binary file if `data.SAVE == True`. The results are plotted if `data.SHOWTRACES == True`.

 The distribution of the synapses is shown if `data.SHOWSYNS == True`

 The cluster size can be changed in lines 92-93:
 240 synapses are clustered, the others are background. Equal-sized clustered are placed randomly on the denritic tree. Number of clusters varies between 4 (60 synapse per cluster) and 240 (1 synapse / cluster)
 
 clust2: 240 = 4x60 = 8x30 = 12x20 = 24x10 = 48x5 = 120x2 = 240x1
 
 Inputs are arranged according to the location of their place fields (CA1) or their orientation preference and then according to their phase preference and tuning width (L2/3). This way consecutive synapses have correlated activity, if connectivity is not randomised.

