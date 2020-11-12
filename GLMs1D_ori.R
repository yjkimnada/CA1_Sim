## A script to simulate population activity ni the visual cortex
## We simulate the response of a population of cells to drifting oriented grating
## cells will have different and orientation and phase preference
## Stimulus: 8 orientations, each shown for 2.5 seconds, grating speed is 2 cycles / sec
## + an additional gamma (40Hz) modulation of the firing rate
## tuning: cells have variable orientation and phase tuning width
## they - are all tuned to the same gamma phase with the same precision
## inhibitory neurons show no orientation or phase tuning, but have a similar tuning to gamma


## distribution of cell-parameters should be matched to V1ori_tuning.key - Neill08.pdf and Durand16.pdf
# 
# monkey - Ecker et al., 2010, Fig S9 - Half-width at half maximum, Orientation tuning index Direction selectivity index, Baseline firing rate, Peak firing rate.

# check orientation
# t > Tmax
# setwd("/data/bbu/V1_Data/")
source('./Functions/gen_Ecells.R', chdir = TRUE)
# source('./gen_Icells.R', chdir = TRUE)
# library('circular')
# require('plotrix')
# library('png')
# library('viridis')
# PAR1 <- 1
# PAR2 <- 1

Ensyn <- 1920; Insyn=192; Ntrial <- 16; Erate <- 5; Irate=40; dt <- 0.001; Tmax=24; Noris=16
Ntrial <- 16
Ntrial <- 2 # comment this out to get the full set of 16 trials
ENs <- c(ori=16, sdOri=5, phase=6, sdPhase=4)
# ENs <- c(ori=16, sdOri=1, phase=6, sdPhase=1)
INs <- c(ori=16, sdOri=1, phase=6, sdPhase=2)
# ENs <- c(ori=1, sdOri=5, phase=6, sdPhase=4); Ntrial <- 2; INs <- c(ori=1, sdOri=1, phase=6, sdPhase=2); Tmax <- 24
types <- c("balanced", "random_NR")
itype <- 2 # 1-3
rseed <- 10
graphics <- F

type <- types[itype]

rand.tuning <- F; rand.rates <- F
if (type == "random_NR") {
	rand.tuning <- T
	rand.rates <- T
}

cat("excitatory simulations started... \n")
sp.E <- gen.V1Ecells(Ns=ENs, Ntrial=Ntrial, rate.ave=Erate, seed=rseed, Noris=Noris, Tmax=Tmax, graphics=graphics, adapt=T, rand.tuning=rand.tuning, rand.rates=rand.rates)
	
cat("inhibitory simulations started... \n")
sp.I <- gen.V1Ecells(Ns=INs, Ntrial= Ntrial, rate.ave=Irate, seed=100 + rseed, graphics=graphics, adapt=T, rand.tuning=rand.tuning, rand.rates=rand.rates, gamma=sp.E$gamma, phases=sp.E$phases, oris=sp.E$oris)
	
# outdir <- paste("~/Projects/KOKI/Synchrony/stims/", type, sep="")
outdir <- paste("./Data/V1/", type, sep="")
dir.create(file.path(outdir), showWarnings = TRUE)
duration <- round(length(sp.E$gamma) * dt)

for (rep in 1:Ntrial){

	Espikes <- sp.E$spt[[rep]] # first column: cell, second column: time in ms
	Espikes[,1] <- Espikes[,1] - 1 # Python index starts from 0
	Espikes[,2] <- round((Espikes[,2] - sp.E$t[rep, 1]) * 1000, 2)
	fname  <- paste (outdir, '/Espikes_d', duration, '_Ne', Ensyn, '_Re', Erate, '_rseed', rseed, '_rep', rep-1, '.dat', sep="")
	write.table(Espikes, file=fname, row.names=FALSE, col.names=FALSE)

	Ispikes <- sp.I$spt[[rep]] # first column: cell, second column: time in ms
	Ispikes[,1] <- Ispikes[,1] - 1
	Ispikes[,2] <- round((Ispikes[,2] - sp.E$t[rep, 1]) * 1000, 2)
	fname  <- paste (outdir, '/Ispikes_d', duration, '_Ni', Insyn, '_Ri', Irate, '_rseed', rseed, '_rep', rep-1, '.dat', sep="")
	write.table(Ispikes, file=fname, row.names=FALSE, col.names=FALSE)
 	
 }
