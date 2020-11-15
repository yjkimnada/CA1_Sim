############################################################
## A script to analyse the dendritic membrane potential responses of the CA1 cell dorung SPWs
## data is recorded at 5000Hz
## Vm in 4 dendritic branches and 4 dendnritic spines, for T=10s and N=16 repetitions
############################################################

source('./functions/scalebar.R')
source('./functions/StatsEval_functions.R', chdir = TRUE)

localdir <- getwd()
# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../CA1/clust3/replay/'
setwd(datadir)

library(viridis)
library(colormap)

stim_type <- 'balanced'
act_type <- 'Dactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)
spinetype <- 'spines'
clust.names <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')

nsim <- length(clust.names)
prop_dspikes <- array(NA, dim=c(2, nsim, 10), dimnames=list(c("control", "LTP"), clust.names, paste("rep", 1:10)))

nseed <- 10
ntrial <- 16
nstart <- 20

Ne <- 2240
Ni <- 200
Re <- 9
Ri <- 30
Tmax <- 0.3
graphics <- FALSE


i.w <- 2; i.clust <- 3; stimseed <- 4
ii.10_11 <- c(1501 * 9 + 500:1000, 1501 * 10 + 500:1000)

for (i.w  in c(1,2)){
	wA <- i.w; wN <- i.w
	clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, spinetype, sep='')
	for (i.clust in 1:7){
		clustSize <- clustSizes[i.clust]
		cat("\n", clustSize, "started \n")
			
		for (stimseed in 1:nseed){
			cat(stimseed, "started \n")
		
			infile <- paste('./', clustSize, '/vDdata_T', Tmax, '_Ne', Ne, '_gA', 0.6, '_tauA', 1, '_gN', 0.8, '_Ni', Ni, '_gG', 0.7, '_gB', 1.2, '_Er', Re, '_Ir', Ri, '_', stim_type, '_rep', ntrial, '_stimseed', stimseed, '.bin', sep="")
			con <- file(infile, "rb")
			dim <- readBin(con, "integer", 2)
			dresp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
			close(con)

			rD <- array(dresp, dim=c(8, 30020, 16))
			# hist(rD[1:4,,], freq=F, col=grey(0.5))
			# hist(rD[1:4, ii.10_11, ], freq=F, col=rgb(1, .5,.5, 0.5), add=T)
			prop_dspikes[i.w, i.clust, stimseed] <- sum(rD[1:4, ii.10_11, ] > -25) / length(rD[1:4, ii.10_11, ])

		}
	}
}

clust_cols <- colormap(colormaps$jet, nshades=30)[c(10, 14, 17, 20, 22, 24, 26)]
LTPcols <- clust_cols
for ( i in 1:7){
	hsvcol <- rgb2hsv(col2rgb(clust_cols)[,i])
	LTP_hsvcol <- hsvcol
	LTP_hsvcol[2] <- LTP_hsvcol[2] * 0.7
	LTP_hsvcol[3] <- LTP_hsvcol[3] * 0.7
	LTPcols[i] <- hsv(LTP_hsvcol[1], LTP_hsvcol[2], LTP_hsvcol[3])
}

nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)


par(mar=c(4,4,1,1))
pplot(data=prop_dspikes[1,,], rep(0, 7), nSynPerClust, cols=clust_cols, pchs=22, xx=1:7*3, ylim=c(0, 1), add=F, greylevel=0.9, xlim=c(2, 23), ylab='dspike ratio')
pplot(data=prop_dspikes[2,,], rep(0, 7), nSynPerClust, cols=LTPcols, pchs=22, xx=1:7*3+1, add=T, greylevel=0.75)
mtext('cluster size', 1, line=2)

dsp_data <- rbind(prop_dspikes[1,,], prop_dspikes[2,,])
rownames(dsp_data) <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60', 'LTP_240_cluster_of_1', 'LTP_120_cluster_of_2', 'LTP_48_cluster_of_5', 'LTP_24_cluster_of_10', 'LTP_12_cluster_of_20', 'LTP_8_cluster_of_30', 'LTP_4_cluster_of_60')
write.table(format(round(dsp_data, 6), scientific=F), file='Fig4g_data.txt')
