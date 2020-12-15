####################################################################################
## A script to analyse the effect of clustering on the postsynaptic Vm response amplitude during SPW
####################################################################################

source('./functions/StatsEval_functions.R', chdir = TRUE)

source('./functions/scalebar.R')
localdir <- getwd()
# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../CA1/clust3/replay/'

library(viridis)
library(colormap)

stim_types <- rep('balanced', 7)
act_type <- 'Dactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)
wA <- 1; wN <- 1
clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')
act_types <- clustSizes
typenames <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')
n.sim <- length(typenames)
graphics <- F


f.in.name <- 'datasets/replay_variability_spines.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data <- eval_replay_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data, file=f.in.name)
}

sds <- stats_data$sds
nSps <- stats_data$nSps
meanresps <- stats_data$meanresps

####################################################################################
wA <- 2; wN <- 2
clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')
act_types <- clustSizes

f.in.name <- 'datasets/replay_variability_spines_LTP.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data_LTP <- eval_replay_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data_LTP, file=f.in.name)
}

sds_LTP <- stats_data_LTP$sds
nSps_LTP <- stats_data_LTP$nSps
meanresps_LTP <- stats_data_LTP$meanresps

####################################################################################
## plots
####################################################################################
clust.sizes <- nSynPerClust
clust_cols <- colormap(colormaps$jet, nshades=30)[c(10, 14, 17, 20, 22, 24, 26)]
LTPcols <- clust_cols
for ( i in 1:7){
	hsvcol <- rgb2hsv(col2rgb(clust_cols)[,i])
	LTP_hsvcol <- hsvcol
	LTP_hsvcol[2] <- LTP_hsvcol[2] * 0.7
	LTP_hsvcol[3] <- LTP_hsvcol[3] * 0.7
	LTPcols[i] <- hsv(LTP_hsvcol[1], LTP_hsvcol[2], LTP_hsvcol[3])
}

########################################################################################
## Fig 4c - spike counts
########################################################################################
spikes <- apply(nSps, c(1,4), mean)
spikes.LTP <- apply(nSps_LTP, c(1,4), mean)

par(mar=c(4,5,1,1))
matplot(t(spikes.LTP), t='l', col=LTPcols, lty=3, axes=F, lwd=3, xlab='trajectory start (cm)', ylab='number of spikes', ylim=c(min(spikes), max(spikes.LTP)));
matplot(t(spikes), t='l', col=clust_cols, lty=1, add=T, lwd=3);
axis(1, seq(1, 20), seq(0, 190, by=10))
axis(2, las=2)

legend('topleft', legend= typenames, lty=1, lwd=3, col=clust_cols, bty='n')
legend('topright', legend= typenames, lty=1, lwd=3, col=LTPcols, bty='n')
# dev.off()

### saving the somatic Vm response amplitudes
ampl <- apply(meanresps[,,150:200,10:11], c(1,2), mean)# - apply(meanresps[,,150:200,16:20], c(1,2), mean)
boxplot(t(ampl), col=clust_cols)
ampl.LTP <- apply(meanresps_LTP[,,150:200,10:11], c(1,2), mean)# - apply(meanresps_LTP[,,150:200,16:20], c(1,2), mean)
boxplot(t(ampl.LTP), col=LTPcols)

save(ampl, file='ampl_SPW.RData')
save(ampl.LTP, file='ampl_SPW_LTP.RData')

########################################################################################
## Fig 4e - SPW amplitude
########################################################################################

baseline <- apply(sds[6,,], 1, mean)
baselineLTP <- apply(sds_LTP[6,,], 1, mean)

par(mar=c(4,4,1,1))
pplot(data=sds[7,,], baseline, nSynPerClust, cols=clust_cols, pchs=22, xx=1:7*3, ylim=c(6, 9), add=F, greylevel=0.9, xlim=c(2, 23), ylab='SPW amplitude (mV)')
pplot(data=sds_LTP[7,,], baselineLTP, nSynPerClust, cols=LTPcols, pchs=22, xx=1:7*3+1, add=T, greylevel=0.75)
mtext('cluster size', 1, line=2)


amp_data <- rbind(sds[7,,], sds_LTP[7,,])
rownames(amp_data) <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60', 'LTP_240_cluster_of_1', 'LTP_120_cluster_of_2', 'LTP_48_cluster_of_5', 'LTP_24_cluster_of_10', 'LTP_12_cluster_of_20', 'LTP_8_cluster_of_30', 'LTP_4_cluster_of_60')

write.table(format(round(amp_data, 4), scientific=F), file='Fig4e_data.txt')

