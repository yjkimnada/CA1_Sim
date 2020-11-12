## Simulate activity during slow wave sleep
##
## we generate 300 ms epochs, 3x100ms: before, during and after
## To do this, we will simulate 60 cm long runs in periods of 100 ms. This equals 6 m/s speed.
##
## The rate: total Nsyn x Nsp x prel = 1600  release events per SPW; 
## 		where nsysn=20000 excitatory synapses, Nsp=0.4 spikes/input and prel=0.2 release probability
## 		this comes to about 8 Hz firing rate at the end
## But the refractoryness decreases the real rate compared to the rate parameter, so we need to increase the rate parameter a bit.
## This script now only works for the Ensyn = 2000, and it takes a few minutes to generate the data.

## First, we simulate trajectory replay in both the clustered and the non-clustered cells

## 1. generate data with different tuning properties
source('./Functions/gen_cells.R', chdir=T) # load the required functions - essentially simulating GLMs in R

Ensyn <- 2000 # number of excitatory inputs to simulate
	# for the clustered cells we separate background activity and replay into different units.
	# the actual number of simulated cells will be Ensyn + Eclustered
Nrep <- 16
Nrep <- 2 ## comment this out to generate all 16 repretitions

Ntrial <- 5*Nrep
Erate <- 9 # rate will be on average 8 Hz - we need to set the parameter much higher because adaptation will decrease it.
iclust <- 881:1120 # the clustered cells have a lower background activity
dt <- 0.001
Erate.bg <- 0.1
Nbasket <- 80 # number of basket cells
Ndend <- 120 # number of dendrite targeting neurons
graphics <- T
types <- c('balanced')
type <- types[1]
sample.template=F; randfields=F; randrates=F; fields.middle=F; same.runs=F
rseed <- 1

####################################
## excitatory cell firing - during SPW
phi.template <- 	c(0,0.4,1.1,0.4)
sp.E <- gen.replay(N.cells=Ensyn, L =200, Ntrial= Ntrial, dt=dt, N.x.basis=40, N.phi.basis=4, rate.ave=Erate, rate.bg=Erate.bg, seed=rseed, mu.v=6, vx=NULL, phi=NULL, graphics=graphics, plot.phase.precess=F, adapt=T, sample.template=sample.template, randfields=randfields, randrates=randrates, fields.middle=fields.middle, var.trials=1/100000, same.runs=same.runs, phi.template=phi.template, i.clust=iclust)

###########################################
## inhibitory cells - perisomatic: 10-30-10 Hz
# during SPWs, high rate, ripple
sp.I_basket <- gen.thetacells(Nbasket, L=200, Ntrial=Ntrial, dt=1/1000, N.phi.basis=4, rate.ave=45, vx=list(x=sp.E$dists, v=sp.E$speed.1D), phi=sp.E$lfp[[1]], seed=23, same.cells=T, phi.template=phi.template, graphics=F)

sp.I_dend <- gen.thetacells(Ndend, L=200, Ntrial=Ntrial, dt=1/1000, N.phi.basis=4, rate.ave=20, vx=list(x=sp.E$dists, v=sp.E$speed.1D), phi=sp.E$lfp[[1]], seed=23, same.cells=T, phi.template=phi.template, graphics=F)


##############################################
# between SPWs
outdir <- "./Data/replay"
dir.create(file.path(outdir), showWarnings = TRUE)
duration <- round(length(sp.E$dists[[1]]) * dt)

tstarts <- seq(0, by=50/3, length=20) * dt # s
dstarts <- seq(0, by=10, length=20)  # cm
iseed <- 1000
Insyn <- Nbasket + Ndend

par(mfcol=c(2,1))
for (ii in 1:20){ # 20 different start positions
	t1 <- tstarts[ii] # start time of the replay sequence
	t2 <- (t1 + 0.1)
	d1 <- dstarts[ii] # start position of the replay
	## we chose the first 16 simulations for ii in c(1, 6, 11, 16)
	kk <- ii %% 5
	for (i.rep in 1: Nrep){ # 16 repetitions
		ind.sim <- kk * Nrep + i.rep
		spE1 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Ensyn, N.c.basis=6, rate.ave=0.8, w.self=NULL, rseed=iseed)$sp
		i.clspikes <- which(spE1[,2] %in% iclust) # distribute the spikes between clustered and background cells
		i.clsp  <- i.clspikes[!!rbinom(length(i.clspikes), 1, 0.9)] # 90% goes to background
		spE1[i.clsp,2] <- spE1[i.clsp,2] + 1120
		
		spE2 <- sp.E$spt[[ind.sim]]
		spE2[,1] <- spE2[,1] - 2*(ind.sim-1)
		spE2[,1] <- round(spE2[,1], 3)
		index.SPW <- which((spE2[,1] > t1) & (spE2[,1] <= t2))
		spE2 <- spE2[index.SPW,]
		spE2[,1] <- spE2[,1] - t1 + 0.1

		nrow(spE2)

		spE3 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Ensyn, N.c.basis=6, rate.ave=0.8, w.self=NULL, rseed=iseed+1)$sp
		i.clspikes <- which(spE3[,2] %in% iclust)
		i.clsp  <- i.clspikes[!!rbinom(length(i.clspikes), 1, 0.9)]
		spE3[i.clsp,2] <- spE3[i.clsp,2] + 1120
		spE3[,1] <- spE3[,1] + 0.2

		Espikes <- rbind(spE1, spE2, spE3)
		Espikes <- Espikes[,c(2,1)]
		Espikes[,2] <- round(Espikes[,2] * 1000,1)
		Espikes[,1] <- round(Espikes[,1] - 1)
		# plot(Espikes[,2], Espikes[,1], pch=16, cex=0.2)

		fname  <- paste (outdir, '/Espikes_d03_Ne', Ensyn, '_Re', Erate, '_dstart', d1, '_rep', i.rep-1, '.dat', sep="")
		write.table(Espikes, file=fname, row.names=FALSE, col.names=FALSE)

		### inhibitory cells

		spIB1 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Nbasket, N.c.basis=6, rate.ave=10.5, w.self=NULL, rseed=iseed+2)$sp
		spID1 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Ndend, N.c.basis=6, rate.ave=5, w.self=NULL, rseed=iseed+4)$sp
		spID1[,2] <- spID1[,2] + Nbasket
		
		spIB2 <- sp.I_basket$spt[[ind.sim]]
		spIB2[,1] <- spIB2[,1] - 2*(ind.sim-1)
		spIB2[,1] <- round(spIB2[,1], 3)
		index.SPW <- which((spIB2[,1] > t1) & (spIB2[,1] <= t2))
		spIB2 <- spIB2[index.SPW,]
		spIB2[,1] <- spIB2[,1] - t1 + 0.1
		
		spID2 <- sp.I_dend$spt[[ind.sim]]
		spID2[,1] <- spID2[,1] - 2*(ind.sim-1)
		spID2[,1] <- round(spID2[,1], 3)
		index.SPW <- which((spID2[,1] > t1) & (spID2[,1] <= t2))
		spID2 <- spID2[index.SPW,]
		spID2[,1] <- spID2[,1] - t1 + 0.1
		spID2[,2] <- spID2[,2] + Nbasket


		spIB3 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Nbasket, N.c.basis=6, rate.ave=10.5, w.self=NULL, rseed=iseed+3)$sp
		spIB3[,1] <- spIB3[,1] + 0.2
		
		spID3 <- gen.adaptcells(Tmax =0.1, dt=0.001, graphics=F, N.cells=Ndend, N.c.basis=6, rate.ave=5, w.self=NULL, rseed=iseed+5)$sp
		spID3[,1] <- spID3[,1] + 0.2
		spID3[,2] <- spID3[,2] + Nbasket


		Ispikes <- rbind(spID1, spID2, spID3, spIB1, spIB2, spIB3)
		tsort <- sort(Ispikes[,1], index=T)$ix
		Ispikes <- Ispikes[tsort,]

		Ispikes <- Ispikes[,c(2,1)]
		Ispikes[,2] <- round(Ispikes[,2] * 1000,1)
		Ispikes[,1] <- round(Ispikes[,1] - 1)

		# plot(Ispikes[,2], Ispikes[,1], pch=16, cex=0.2)

		fname  <- paste (outdir, '/Ispikes_d03_Ni', Insyn, '_Ri30_dstart', d1, '_rep', i.rep-1, '.dat', sep="")
		write.table(Ispikes, file=fname, row.names=FALSE, col.names=FALSE)

		iseed <- iseed + 10
		cat(ii, i.rep, '\n')
	}
	
}

