########################################################################
## A script to analyse the effect of the presence of strong dendritic branches and synaptic cluatering on the somatic response
########################################################################

library(viridis)
library(colormap)
library(gplots)

source('./functions/StatsEval_functions.R', chdir=T)

localdir <- getwd()
# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../CA1/clust3/NaStrong/'
setwd(localdir)

clust.names <- c(1, 2, 5, 10, 20, 30, 60)


stim_types <- rep('balanced', 7)
act_type <- 'Bactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)
wA <- 1; wN <- 1
act_types <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')
typenames <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')
n.sim <- length(typenames)
graphics <- F


f.in.name <- 'datasets/replay_variability_trial1.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data <- eval_strong_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data, file=f.in.name)
}

ampl_base <- stats_data$ampl_base
ampl_clust <- stats_data$ampl_clust
nSps <- stats_data$nSps
meanresps <- stats_data$meanresps

####################################################################################
wA <- 2; wN <- 2
act_types <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')

f.in.name <- 'datasets/replay_variability_trial1_LTP.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data_LTP <- eval_strong_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data_LTP, file=f.in.name)
}

ampl_base_LTP <- stats_data_LTP$ampl_base
ampl_clust_LTP <- stats_data_LTP$ampl_clust
nSps_LTP <- stats_data_LTP$nSps
meanresps_LTP <- stats_data_LTP$meanresps

####################################################################################
wA <- 1; wN <- 1
act_types <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines_hotspot', sep='')

f.in.name <- 'datasets/replay_variability_trial1_hotspot.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data_hotspot <- eval_strong_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data_hotspot, file=f.in.name)
}

ampl_base_hotspot <- stats_data_hotspot$ampl_base
ampl_clust_hotspot <- stats_data_hotspot$ampl_clust
nSps_hotspot <- stats_data_hotspot$nSps
meanresps_hotspot <- stats_data_hotspot$meanresps

####################################################################################
wA <- 2; wN <- 2
act_types <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines_hotspot', sep='')

f.in.name <- 'datasets/replay_variability_trial1_LTP_hotspot.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data_hotspot_LTP <- eval_strong_stats(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=5000)
	setwd(localdir)
	save(stats_data_hotspot_LTP, file=f.in.name)
}

ampl_base_hotspot_LTP <- stats_data_hotspot_LTP$ampl_base
ampl_clust_hotspot_LTP <- stats_data_hotspot_LTP$ampl_clust
nSps_hotspot_LTP <- stats_data_hotspot_LTP$nSps
meanresps_hotspot_LTP <- stats_data_hotspot_LTP$meanresps


####################################################################################
## plots
####################################################################################
setwd('~/Projects/KOKI/Synchrony/CA1/clust3/NaStrong')
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
## Fig 5c - spike counts
########################################################################################
i.start <- 9
par(mfcol=c(2,1))
par(mar=c(1,4,4,1))

mm <- apply(nSps[,,i.start], 1, mean)
ss <- apply(nSps[,,i.start], 1, sd)
plotCI(0:6*2 + 1/2, mm, ss, cex=2, gap=0, pch=21, pt.bg=clust_cols, axes=F, xlab='', ylab='spike count', xlim=c(-1, 14), ylim=c(0, 5), t='o')

mm <- apply(nSps_hotspot[,,i.start], 1, mean)
ss <- apply(nSps_hotspot[,,i.start], 1, sd)
plotCI(0:6*2 + 1, mm, ss, cex=2, gap=0, col=clust_cols, add=T)
axis(2, las=2)
white <- grey(1)
legend('topleft', leg=c('control', 'LTP', 'hotspot', 'hotspot+LTP'), pch=21, pt.bg=c(clust_cols[5], LTPcols[5], white, white), col=c(1,1,clust_cols[5], LTPcols[5]), bty='n')


par(mar=c(4,4,1,1))
mm <- apply(nSps_LTP[,,i.start], 1, mean)
ss <- apply(nSps_LTP[,,i.start], 1, sd)
plotCI(0:6*2 + 1/2, mm, ss, cex=2, gap=0, pch=21, pt.bg=LTPcols, axes=F, t='o', xlab='cluster size', ylab='spike count', xlim=c(-1, 14), ylim=c(0, 5))

mm <- apply(nSps_hotspot_LTP[,,i.start], 1, mean)
ss <- apply(nSps_hotspot_LTP[,,i.start], 1, sd)
plotCI(0:6*2 + 1, mm, ss, cex=2, gap=0, col=LTPcols, add=T)

axis(1, 0:6*2+3/10, clust.names, tick=F)
axis(2, las=2)
# dev.off()


spcount_data <- rbind(nSps[,,9], nSps_hotspot[,,9], nSps_LTP[,,9], nSps_hotspot_LTP[,,9])
rnames <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60')
rownames(spcount_data) <- paste(rep(c('', 'hotspot_', 'LTP_', 'hotspot_LTP_'), e=7), rnames, sep='')
write.table(spcount_data, file='Fig5d_data.txt')


########################################################################################
## Fig 5d - size of SPW amplitude
########################################################################################

par(mfcol=c(2,1))
par(mar=c(1,4,4,1))

mm <- apply(ampl_clust[,,i.start]-ampl_base[,,i.start], 1, mean)
ss <- apply(ampl_clust[,,i.start]-ampl_base[,,i.start], 1, sd)
plotCI(0:6*2 + 1/2, mm, ss, cex=2, gap=0, pch=21, pt.bg=clust_cols, axes=F, xlab='', ylab='SPW amplitude (mV)', xlim=c(-1, 14), ylim=c(5, 9), t='o')

mm <- apply(ampl_clust_hotspot[,,i.start]-ampl_base_hotspot[,,i.start], 1, mean)
ss <- apply(ampl_clust_hotspot[,,i.start]-ampl_base_hotspot[,,i.start], 1, sd)
plotCI(0:6*2 + 1, mm, ss, cex=2, gap=0, col=clust_cols, add=T)
axis(2, las=2)
white <- grey(1)
legend('topleft', leg=c('control', 'LTP', 'hotspot', 'hotspot+LTP'), pch=21, pt.bg=c(clust_cols[5], LTPcols[5], white, white), col=c(1,1,clust_cols[5], LTPcols[5]), bty='n')


par(mar=c(4,4,1,1))
mm <- apply(ampl_clust_LTP[,,i.start]-ampl_base_LTP[,,i.start], 1, mean)
ss <- apply(ampl_clust_LTP[,,i.start]-ampl_base_LTP[,,i.start], 1, sd)
plotCI(0:6*2 + 1/2, mm, ss, cex=2, gap=0, pch=21, pt.bg=LTPcols, axes=F, t='o', xlab='cluster size', ylab='SPW amplitude (mV)', xlim=c(-1, 14), ylim=c(5, 9))

mm <- apply(ampl_clust_hotspot_LTP[,,i.start]-ampl_base_hotspot_LTP[,,i.start], 1, mean)
ss <- apply(ampl_clust_hotspot_LTP[,,i.start]-ampl_base_hotspot_LTP[,,i.start], 1, sd)
plotCI(0:6*2 + 1, mm, ss, cex=2, gap=0, col=LTPcols, add=T)

axis(1, 0:6*2+3/10, clust.names, tick=F)
axis(2, las=2)

SPWampl_data <- rbind(ampl_clust[,,9]-ampl_base[,,9], ampl_clust_hotspot[,,9]-ampl_base_hotspot[,,9], ampl_clust_LTP[,,9]-ampl_base_LTP[,,9], ampl_clust_hotspot_LTP[,,9]-ampl_base_hotspot_LTP[,,9])
rnames <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60')
write.table(format(round(SPWampl_data, 4), scientific=F), file='Fig5c_data.txt')

