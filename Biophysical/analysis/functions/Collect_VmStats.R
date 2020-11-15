##########################################
## A script to analyse the effect of various factors on the postsynaptic Vm
##########################################

source('./StatsEval_functions.R', chdir = TRUE)

## data: binary file, containing the somatic membrane potential of the biophysical model
## scripts for simulating the biophysical model can be found here: https://bitbucket.org/bbu20/clustering/
## data format: matrix, with rows being the individual trials, recorded at 1000 Hz
localdir <- getwd()

# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../../CA1/global_regular_lrand/place/'
library(viridis)
library(colormap)
library(gplots)

stim_types <- c('random_NR', 'balanced')
act_types <- c('Dactive', 'Dactive')
spinetype <- 'spines'
typenames <- c('act_rand', 'act_bal')
n.sim <- length(stim_types)

ntrial <- 16
Tmax <- 10
gN <- '0.8'

graphics <- F

fname <- '../datasets/response_variability_ttt_globalReg_localRand_spines.Rdata'
if (file.exists(fname)) load(fname) else {
	setwd(datadir)
	stats_data <- eval_Vm_stats(n.sim, stim_types, act_types, typenames, graphics=0)
	setwd(localdir)
	save(stats_data, file=fname)
}

sds <- stats_data$sds
meanresps <- stats_data$meanresps
ttt_variance <- stats_data$ttt_variance
ttt_variance_peak <- stats_data$ttt_variance_peak

############################################################
cols <- c(rgb(255, 208, 220, max=255), rgb(181, 181, 181, max=255))
pchs <- c(22, 22, 23, 23)

plot_data_SEM(filename=NULL, sds, ttt_variance, ttt_variance_peak, typenames, col=cols, pch=pchs)~/Projects/KOKI/Synchrony/CA1/
