##########################################
## A script to analyse the NMDA current during clustering
##########################################
library(viridis)
library(gplots)
library(colormap)

localdir <- getwd()


stim_type <- 'balanced'
act_type <- 'Dactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)

clust.names <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')

spinetype <- 'spines'
n.reps <- 10
nsim <- length(clust.names)

ww <- 2; i.clust <- 1; stimseed <- 1
INs <- array(NA, dim=c(2, 7, 40, 20))


# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
place <- T
if (place) {
	datadir <- '../CA1/clust2/place/' 
} else {
	datadir <- '../CA1/clust3/replay/'
}
setwd(datadir)

for (ww in 1:2){
	wA <- ww; wN <- ww
	clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, spinetype, sep='')
			
	for (i.clust in 1:nsim){
		clustSize <- clustSizes[i.clust]
		
		k <- 0
		for (stimseed in 1:n.reps){
			# cat(stimseed, "started \n")
			infile <- paste('./', clustSize, '/mean_IN_', stimseed, '.bin', sep='')
			con <- file(infile, "rb")
			dim <- readBin(con, "integer", 2)
			resp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
			close(con)
			INs[ww, i.clust,1:4 + k,] <- resp
			k <- k + 4
		}
	}
}

setwd(localdir)

if (place) {
	save(INs, file='./datasets/INs_place.RData')
} else {
	save(INs, file='./datasets/INs_SPW.RData')
}

################################################################

load(file='datasets/INs_SPW.RData')
Inmda_r <- INs * 1000

load(file='datasets/INs_place.RData')
Inmda_p <- INs * 1000


## plot the NMDA current at the different parts of the trajectory
## symbols mark the analysed data
clust_cols <- colormap(colormaps$jet, nshades=30)[c(10, 14, 17, 20, 22, 24, 26)]
LTPcols <- clust_cols
for ( i in 1:7){
	hsvcol <- rgb2hsv(col2rgb(clust_cols)[,i])
	LTP_hsvcol <- hsvcol
	LTP_hsvcol[2] <- LTP_hsvcol[2] * 0.7
	LTP_hsvcol[3] <- LTP_hsvcol[3] * 0.7
	LTPcols[i] <- hsv(LTP_hsvcol[1], LTP_hsvcol[2], LTP_hsvcol[3])
}

mIN_p <- apply(Inmda_p, c(1,2,4), mean)
mIN_r <- apply(Inmda_r, c(1,2,4), mean)

matplot(t(mIN_p[2,,]), t='l', lty=1, col=LTPcols, ylim=c(-6, 6), xlab='trajectory start', ylab='NMDA current (pA)')
matplot(t(mIN_p[1,,]), t='l', lty=1, col=clust_cols, ylim=range(mIN_p), add=T)
matplot(c(9, 10), t(mIN_p[2,,c(9, 10)]), pch=21, bg=LTPcols, add=T, col=1)
matplot(c(9, 10), t(mIN_p[1,,c(9, 10)]), pch=21, bg=clust_cols, add=T, col=1)


matplot(-1 * t(mIN_r[2,,]), t='l', lty=2, col=LTPcols, ylim=range(mIN_r), add=T, lwd=2)
matplot(-1 * t(mIN_r[1,,]), t='l', lty=2, col=clust_cols, ylim=range(mIN_r), add=T, lwd=2)
matplot(c(10, 11), -1 * t(mIN_r[2,,c(10, 11)]), pch=21, bg=LTPcols, add=T, col=1)
matplot(c(10, 11), -1*t(mIN_r[1,,c(10, 11)]), pch=21, bg=clust_cols, add=T, col=1)

### the extra current due to clustering is the difference between the currentin the clustered versus non-clustered configuration
IN_p7 <- apply(Inmda_p[,,,9:10], c(1,2,3), mean)
IN_p6 <- IN_p7[,2:7,]
for (j in 1:2){
	for (i in 1:6) IN_p6[j,i,] <- IN_p6[j,i,] - IN_p7[1,1,]
}

IN_r7 <- apply(Inmda_r[,,,10:11], c(1,2,3), mean)
IN_r6 <- IN_r7[,2:7,]
for (j in 1:2){
	for (i in 1:6) IN_r6[j,i,] <- IN_r6[j,i,] - IN_r7[1,1,]
}


col.SPW <- rgb(58, 144, 255, max=255)
col.SPW_LTP <- rgb(203, 63, 219, max=255)
col.theta <- rgb(54, 92, 141, max=255)
col.theta_LTP <- rgb(101, 21, 110, max=255)
fill.theta <- rgb(31, 161, 135, max=255)
fill.theta_LTP <- rgb(212, 72, 66, max=255)

##########################################################################
### Fig 4d bottom
##########################################################################

par(mar=c(4,4,1,1))
m1 <- apply(IN_p7[1,,], 1, mean); s1 <- apply(IN_p7[1,,], 1, sd) / sqrt(40)
m2 <- apply(IN_p7[2,,], 1, mean); s2 <- apply(IN_p7[2,,], 1, sd) / sqrt(40)
plotCI(1:7, m1, s1, pt.bg=fill.theta, col=col.theta, t='o', pch=21, lty=1, ylim=c(-6, 0), axes=F, xlab='cluster size', ylab='I_NMDA, pA', main='', gap=0)
plotCI(1:7, m2, s2, pt.bg=fill.theta_LTP, col=col.theta_LTP, t='o', pch=21, lty=1, add=T, gap=0)
axis(2, las=2); axis(1, 1:7, nSynPerClust)


m1 <- apply(IN_r7[1,,], 1, mean); s1 <- apply(IN_r7[1,,], 1, sd) / sqrt(40)
m2 <- apply(IN_r7[2,,], 1, mean); s2 <- apply(IN_r7[2,,], 1, sd) / sqrt(40)
plotCI(1:7, m1, s1, col=col.SPW, t='o', pch=21, pt.bg=grey(1), lty=1, add=T, gap=0)
plotCI(1:7, m2, s2, col=col.SPW_LTP, t='o', pch=21, pt.bg=grey(1), lty=1, add=T, gap=0)
axis(2, las=2); axis(1, 1:7, nSynPerClust)



iNMDA_data <- rbind(IN_p7[1,,], IN_p7[2,,], IN_r7[1,,], IN_r7[2,,])
rnames <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60')
rownames(iNMDA_data) <- paste(rep(c('iN_place', 'iN_place_LTP', 'iN_SPW', 'iN_SPW_LTP'), e=7), rnames, sep='')
write.table(format(round(iNMDA_data, 3), scientific=F), file='Fig4d_top_data.txt')

##########################################################################
### Fig 4d top
##########################################################################

## mean depolarisation of the 10-11th sVm during the SPW, 
ampl_r <- array(NA, dim=c(2, 7, 40))
load('datasets/ampl_SPW_LTP.RData')
ampl_r[2,,] <- matrix(rep(t(ampl.LTP), each=4), 7, byrow=T)
load('datasets/ampl_SPW.RData')
ampl_r[1,,] <- matrix(rep(t(ampl), each=4), 7, byrow=T)

## mean depolarisation of the 4-5s of the sVm during running, 9-10 bins
ampl_p <- array(NA, dim=c(2, 7, 40))
load('datasets/ampl_place_LTP.RData')
ampl_p[2,,] <- matrix(rep(t(ampl_LTP), each=4), 7, byrow=T)
load('datasets/ampl_place.RData')
ampl_p[1,,] <- matrix(rep(t(ampl), each=4), 7, byrow=T)


ampl_p6 <- ampl_p[,2:7,]
for (j in 1:2){
	for (i in 1:6) ampl_p6[j,i,] <- ampl_p6[j,i,] - ampl_p[1,1,]
}

ampl_r6 <- ampl_r[,2:7,]
for (j in 1:2){
	for (i in 1:6) ampl_r6[j,i,] <- ampl_r6[j,i,] - ampl_r[1,1,]
}


linreg2_r <- lm(as.vector(apply(ampl_r6, c(1,2), mean)) ~ as.vector(apply(IN_r6, c(1,2), mean)))
linreg2_p <- lm(as.vector(apply(ampl_p6, c(1,2), mean)) ~ as.vector(apply(IN_p6, c(1,2), mean)))

mI_r <- c(apply(IN_r6[1,,], 1, mean), apply(IN_r6[2,,], 1, mean))
sI_r <- c(apply(IN_r6[1,,], 1, sd), apply(IN_r6[2,,], 1, sd))# / sqrt(40)
mA_r <- c(apply(ampl_r6[1,,], 1, mean), apply(ampl_r6[2,,], 1, mean))
sA_r <- c(apply(ampl_r6[1,,], 1, sd), apply(ampl_r6[2,,], 1, sd))# / sqrt(40)

mI_p <- c(apply(IN_p6[1,,], 1, mean), apply(IN_p6[2,,], 1, mean))
sI_p <- c(apply(IN_p6[1,,], 1, sd), apply(IN_p6[2,,], 1, sd))# / sqrt(40)
mA_p <- c(apply(ampl_p6[1,,], 1, mean), apply(ampl_p6[2,,], 1, mean))
sA_p <- c(apply(ampl_p6[1,,], 1, sd), apply(ampl_p6[2,,], 1, sd))# / sqrt(40)


################################################################
## Fig4h
################################################################
plotCI(x=mI_r, y=mA_r, uiw=sA_r, gap=0.5, xlab='NMDA current (pA)', ylab='response amplitude (mV)', axes=F, xlim=c(-5, 1), ylim=c(-.5, 6), pch=c(1,1,1,1,1,1,2,2,2,2,2,2), col=c(clust_cols[2:7], LTPcols[2:7]), cex=1)
plotCI(x=mI_r, y=mA_r, uiw=sI_r, err='x', gap=0.5, pch=c(1,1,1,1,1,1,2,2,2,2,2,2), col=c(clust_cols[2:7], LTPcols[2:7]), cex=1, add=T)

points(mI_p, mA_p, pch=c(1,1,1,1,1,1,2,2,2,2,2,2)+20, bg=c(clust_cols[2:7], LTPcols[2:7]), cex=1, add=T)
plotCI(mI_p, mA_p, uiw=sA_p, gap=0.5, pch='', col=c(clust_cols[2:7], LTPcols[2:7]), cex=1, add=T)
plotCI(mI_p, mA_p, uiw=sI_p, err='x', gap=0.5, pch='', col=c(clust_cols[2:7], LTPcols[2:7]), cex=1, add=T)
axis(1); axis(2, las=2)

abline(linreg2_p, col=1)
abline(linreg2_r, lty=3, col=1)
abline(h=0, v=0, col=grey(0.75))

legend('topright', leg= c(paste(nSynPerClust[-1], 'per clust'), 'LTP', 'place', 'replay'), pch=c(rep(21,9), 1), pt.bg=c(clust_cols[-1], LTPcols[6], clust_cols[6], NA), col=c(rep(1, 9), clust_cols[6]),  bty='n')


cor_data <- rbind(IN_p6[1,,], IN_p6[2,,], IN_r6[1,,], IN_r6[2,,], ampl_p6[1,,], ampl_p6[2,,], ampl_r6[1,,], ampl_r6[2,,])
rnames <- c('120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60')
rownames(cor_data) <- paste(rep(c('iN_place', 'iN_place_LTP', 'iN_SPW', 'iN_SPW_LTP', 'ampl_place', 'ampl_place_LTP', 'ampl_SPW', 'ampl_SPW_LTP'), e=6), rnames, sep='')
write.table(format(round(cor_data, 3), scientific=F), file='Fig4h_data.txt')

