### General information

Library for simulating synthetic hippocampal population activity

for more details, see our paper on bioRxive:
**Balazs B Ujfalussy and Judit K Makara:** *Impact of functional synapse clusters on neuronal response selectivity* `https://doi.org/ 10.1101/634220`

author: Balazs B Ujfalussy, 2019 `balazs.ujfalussy@gmail.com`

This code is required to reproduces some of the data in our paper: the in vivo-like stimuli for CA1 and L2/3 cells, which we used to test the cluster sensitivity of the neurons.

### Requirements

To run the code, you need R and a few packages (Matrix, viridis, png) installed. 

The code has been tested with R 3.4.0, installed on a Macbook Pro Mid 2015, macOS Mojave 10.14.5, 

### Details

There are four conditions, they can be generated using different files:
GLMs1D_place.R:	Theta - balanced or random

GLMs1D_replay.R: Sharp Wave - balanced

GLMs1D_ori.R: V1 - balanced or random


### Expected output
Output is saved in subfolders Place, Replay and V1 in the Data folder. Separate output is saved for inhibitory and excitatory cells.

The output format is the following: 

Espikes_D10_Ne2000_Re1_rseed1_rep14.dat:
	excitatory spike times for 10s and for 2000 cells with 1 Hz average firing rate and random seed 1 and repetition 14. 
	First column is cell id, second column is spike time (ms)

Ispikes_D10_Ne200_Re7.5_rseed1_rep12.dat:
	excitatory spike times for 10s and for 2000 cells with 1 Hz average firing rate and random seed 1 and repetition 14. 
	First column is cell id, second column is spike time (ms)

rseed: random seed, controlling the selection of the presynaptic place fields (how many cells active at any given location)
rep: trials with identical presynaptic firing rates but different Poisson spikes
