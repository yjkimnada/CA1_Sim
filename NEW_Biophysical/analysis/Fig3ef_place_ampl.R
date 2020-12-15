#############################################################
## A script to analyse the effect of clustering on the postsynaptic Vm
#############################################################

source('./functions/StatsEval_functions.R', chdir = TRUE)
source('./functions/scalebar.R')
localdir <- getwd()
# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../CA1/clust2/place/'

library(viridis)
library(colormap)
library(gplots)

stim_types <- rep('balanced', 7)
act_type <- 'Dactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)
wA <- 1; wN <- 1
clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')
act_types <- clustSizes
clust.names <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')
n.sim <- length(clust.names)
graphics <- F


f.in.name <- 'datasets/response_variability_spines.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data <- eval_Vm_stats(n.sim, stim_types, act_types, clust.names, graphics=0, ramp=T, sampling_freq=5000)
	setwd(localdir)
	save(stats_data, file=f.in.name)
}

# f.in.name <- 'response_variability_spines_largeG.Rdata'
# if (file.exists(f.in.name)) load(f.in.name) else {
	# stats_data <- eval_Vm_stats(n.sim, stim_types, act_types, clust.names, graphics=0, ramp=T, sampling_freq=5000)
	# save(stats_data, file=f.in.name)
# }

sds <- stats_data$sds
meanresps <- stats_data$meanresps
ttt_variance <- stats_data$ttt_variance
ttt_variance_peak <- stats_data$ttt_variance_peak


###################################

wA <- 2; wN <- 2
clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, 'spines', sep='')
act_types <- clustSizes

f.in.name <- 'datasets/response_variability_LTPspines.Rdata'
if (file.exists(f.in.name)) load(f.in.name) else {
	setwd(datadir)
	stats_data_LTP <- eval_Vm_stats(n.sim, stim_types, act_types, clust.names, graphics=0, ramp=T, sampling_freq=5000)
	setwd(localdir)
	save(stats_data_LTP, file=f.in.name)
}


sds_LTP <- stats_data_LTP$sds
meanresps_LTP <- stats_data_LTP$meanresps
ttt_variance_LTP <- stats_data_LTP$ttt_variance
ttt_variance_peak_LTP <- stats_data_LTP$ttt_variance_peak

###################################

clust_cols <- colormap(colormaps$jet, nshades=30)[c(10, 14, 17, 20, 22, 24, 26)]
LTPcols <- clust_cols
for ( i in 1:7){
	hsvcol <- rgb2hsv(col2rgb(clust_cols)[,i])
	LTP_hsvcol <- hsvcol
	LTP_hsvcol[2] <- LTP_hsvcol[2] * 0.7
	LTP_hsvcol[3] <- LTP_hsvcol[3] * 0.7
	LTPcols[i] <- hsv(LTP_hsvcol[1], LTP_hsvcol[2], LTP_hsvcol[3])
}

pchs <- c(22, 22, 22, 22, 22, 22, 22)
# fname <- 'response_variability_spines.pdf'
# plot_data(filename=NULL, sds, ttt_variance, ttt_variance_peak, clust.names, col=clust_cols, pch=pchs)
# plot_data(filename=fname, sds, ttt_variance, ttt_variance_peak, clust.names, col=clust_cols, pch=pchs)


### plot all the different statistics in the same plot
xx <- seq(0, 18, by=3)
plot_data_SEM_clust(filename=NULL, sds, ttt_variance, ttt_variance_peak, clust.names, col= clust_cols, pch=pchs, xx=xx)

xx <- seq(1, 19, by=3)
plot_data_SEM_clust(filename=NULL, sds_LTP, ttt_variance_LTP, ttt_variance_peak_LTP, clust.names, col= LTPcols, pch=pchs, xx=xx)

###########################################
## Fig 3e
###########################################

par(mfcol=c(1,1))
par(mar=c(2,4,2,1))

data <- sds[1,,]^2
baseline <- apply(ttt_variance, 1, mean) / 16
xlabs <- c(1,2,5,10, 20, 30, 60)
pplot(filename=NULL, data, baseline, xlabs, cols=clust_cols, pchs=22, xx=seq(0, 18, by=3), xlim=c(-1, 20), ylim=c(0, 3), main='TC variance', ylab='tuning curve variance (mV^2)')
lines(seq(0, 18, by=3), rowMeans(data))

data <- sds_LTP[1,,]^2
baseline <- apply(ttt_variance_LTP, 1, mean) / 16
xlabs <- c(1,2,5,10, 20, 30, 60)
pplot(filename=NULL, data, baseline, xlabs, cols=LTPcols, pchs=22, xx=seq(0, 18, by=3)+1, xlim=c(-1, 20), add=T)
abline(h=0, lty=3)
lines(seq(1, 19, by=3), rowMeans(data))


TCvar_data <- rbind(sds[1,,]^2, sds_LTP[1,,]^2)
rownames(TCvar_data) <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60', 'LTP_240_cluster_of_1', 'LTP_120_cluster_of_2', 'LTP_48_cluster_of_5', 'LTP_24_cluster_of_10', 'LTP_12_cluster_of_20', 'LTP_8_cluster_of_30', 'LTP_4_cluster_of_60')
write.table(format(round(TCvar_data, 6), scientific=F), file='Fig3e_data.txt')


###########################################
## Fig 3f
###########################################
par(mfcol=c(1,1))
par(mar=c(2,4,2,1))

data <- sds[6,,]
baseline <- rep(0, 7)
xlabs <- c(1,2,5,10, 20, 30, 60)
pplot(filename=NULL, data, baseline, xlabs, cols=clust_cols, pchs=22, xx=seq(0, 18, by=3), xlim=c(-1, 20), ylim=c(-1, 10), greylevel=0.88, main='TC integral', ylab='response integral (mV x s)')
lines(seq(0, 18, by=3), rowMeans(data))
# dev.off()

data <- sds_LTP[6,,]
baseline <- rep(0, 7)
xlabs <- c(1,2,5,10, 20, 30, 60)
pplot(filename=NULL, data, baseline, xlabs, cols=LTPcols, pchs=22, xx=seq(0, 18, by=3)+1, xlim=c(-1, 20), add=T)
abline(h=0, lty=3)
lines(seq(1, 19, by=3), rowMeans(data))
# dev.off()



respint_data <- rbind(sds[6,,], sds_LTP[6,,])
rownames(respint_data) <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60', 'LTP_240_cluster_of_1', 'LTP_120_cluster_of_2', 'LTP_48_cluster_of_5', 'LTP_24_cluster_of_10', 'LTP_12_cluster_of_20', 'LTP_8_cluster_of_30', 'LTP_4_cluster_of_60')
write.table(format(round(respint_data, 6), scientific=F), file='Fig3f_data.txt')


###########################################
## Fig 3d
###########################################
######################################################
 # look at the mean responses - and save the amplitudes!
 
mresp <- apply(meanresps, c(1,3), mean, na.rm=T)
mresp_LTP <- apply(meanresps_LTP, c(1,3), mean, na.rm=T)

par(mar=c(1,1,1,1))
lty <- c(1,1,1,1, 1,1,1)
matplot(t(mresp), t='l', lwd=c(2), col=clust_cols, axes=F, ylim=c(-68, -62), lty=lty, xlab='', ylab='')
matplot(t(mresp_LTP), t='l', lwd=c(2), col=LTPcols, axes=F, ylim=c(-67.5, -62), lty=lty, xlab='', ylab='', add=T)
legend('topleft', legend=clust.names, lty=lty, lwd=c(2), col=clust_cols, bty='n')
scalebar2(1000, 0.5, '1 s', '0.5 mV', 'topright')

###########################################
## analysing response amplitudes for Fig 4h
###########################################

base <- apply(meanresps[,,1500:2500], c(1,2), mean)
maxV <- apply(meanresps[,,4000:5000], c(1,2), mean)
ampl <- maxV# - base

par(mar=c(4,4,1,1))
matplot(ampl, t='o', pch=16, col=viridis(10)[7], lty=1)
lines(rowMeans(ampl), lwd=3)
save(ampl, file='datasets/ampl_place.RData')


base <- apply(meanresps_LTP[,,1000:2000], c(1,2), mean)
maxV <- apply(meanresps_LTP[,,4000:5000], c(1,2), mean)
ampl_LTP <- maxV# - base

matplot(ampl_LTP, t='o', pch=16, col=viridis(10)[7], lty=1, ylim=c(-68, -61))
matplot(ampl, t='o', pch=16, col=viridis(10)[9], lty=1, add=T)
lines(rowMeans(ampl), lwd=3)
lines(rowMeans(ampl_LTP), lwd=3)
save(ampl_LTP, file='datasets/ampl_place_LTP.RData')

