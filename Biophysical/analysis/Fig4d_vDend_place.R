###################################################
## A script to analyse dendritic membrane potentials during theta / place cells
###################################################

# datadir is the directory where the biophysical model's output was saved 
# it is now pointing to a non-existing directory, which should be created when running the biophysical model with appropriate parameters.
datadir <- '../CA1/clust2/place/'
localdir <- getwd()
setwd(datadir)

library(viridis)
library(colormap)

stim_type <- 'balanced'
act_type <- 'Dactive'
nClust <- c(240, 120, 48, 24, 12, 8, 4)
nSynPerClust <- c(1, 2, 5, 10, 20, 30, 60)

brs <- -100:50
vvv <- -101:50 + 0.5
vv <- -100:49 + 0.5
dvv <- -99:49
JM <- function(v){1 / (1 + exp(-0.071 * v) * (1 / 4.3))}
g_NMDA <- JM(vv) * 0.6
I_NMDA <- g_NMDA * vv
I_NMDA2 <- JM(dvv) * 0.6 * dvv

dI_dv <- diff(JM(brs) * 0.6 * brs)
d2I_dv2 <- diff(diff(JM(vvv) * 0.6 * vvv))

clust.names <- c('240 cluster of 1', '120 cluster of 2', '48 cluster of 5', '24 cluster of 10', '12 cluster of 20', '8 cluster of 30', '4 cluster of 60')
clust_index <- array(NA, dim=c(2, 7, 10, 4), dimnames=list(paste('w', 1:2), clust.names, paste('seed', 1:10), paste('dend', 1:4)))
E_Vd <- array(NA, dim=c(2, 7, 10, 4), dimnames=list(paste('w', 1:2), clust.names, paste('seed', 1:10), paste('dend', 1:4)))
E_I <- array(NA, dim=c(2, 7, 10, 4), dimnames=list(paste('w', 1:2), clust.names, paste('seed', 1:10), paste('dend', 1:4)))
E_I_rel <- array(NA, dim=c(2, 7, 10, 4), dimnames=list(paste('w', 1:2), clust.names, paste('seed', 1:10), paste('dend', 1:4)))

layout(matrix(c(1:7*2, 1:7*2-1, 1:7*2-1+14, 1:7*2+14), 7), heights = c(3,3,3,3,3,3,4))
par(mar=c(0.5,3,1,1))
ww <- 2; i.clust <- 1; stimseed <- 1

for (ww in 1:2){
	wA <- ww; wN <- ww
	spinetype <- 'spines'
	clustSizes <- paste(act_type, '_midMaze_', nClust, 'clusterOf', nSynPerClust, '_wA', wA, '_wN', wN, spinetype, sep='')
		
	nsim <- length(clust.names)
	n.reps <- 10
	
	for (i.clust in 1:nsim){
		clustSize <- clustSizes[i.clust]
		hvDs <- array(NA, dim=c(150, 40, 20)) # 1mV voltage bins - 10 trials and 4 dendrites - 20 spatial bins
		
		k <- 1
		for (stimseed in 1:n.reps){
			# cat(stimseed, "started \n")
			infile <- paste('./', clustSize, '/vDdata_T10_Ne2000_gA0.6_tauA1_gN0.8_Ni200_gG0.7_gB1.2_Er1_Ir7.5_balanced_rep16_stimseed', stimseed, '.bin', sep='')
			con <- file(infile, "rb")
			dim <- readBin(con, "integer", 2)
			resp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
			close(con)
			for (i.branch in 1:4){
				v.branch <- matrix(resp[i.branch,], ncol=16)
				for (i.space in 1:20){
					ii.spbins <- 1:2500 + (i.space-1) * 2500
					hvDs[,k,i.space] <- hist(v.branch[ii.spbins,], br=brs, plot=F)$counts
				}
				k <- k + 1
			}
		}
	
		hvD <- apply(hvDs, c(1,3), sum)
	
		# hist(vDs, br=brs, main='')
		if (i.clust == 1){
			if (ww==1) main <- 'control' else main <- 'LTP'
		} else {
			main <- ''
		}
		if (i.clust == 7){
			par(mar=c(4,3.5,1,1))
		} else {
			par(mar=c(0.5,3.5,1,1))
		}

		image(1:20, vv, t(sqrt(hvD)),col=viridis(40), axes=F, xlab='', ylab='', main=main, ylim=c(-80, 0))
		abline(v=9, col=rgb(1, 0.6, 0))
		if (i.clust == 7){
			axis(2, las=2); axis(1, 0:4*5+0.5, 0:4*50 );
			mtext('dendritic voltage (mV)', 2, 2.5, cex=0.7)
			mtext('position (cm)', 1, 2, cex=0.7)
		}
		# hist_vDs[ww, i.clust,,] <- hvD
		
		vDs <- hvDs[,,9]+hvDs[,,10]
		pvDs <- vDs / colSums(vDs)
		dpvDs <- (-1) * apply(pvDs, 2, diff)
		# Rin <- as.vector(R_ins[ww, i.clust,,])
		E_Is <- round(colSums(pvDs * I_NMDA * ww),3)

		vDs_rel <- hvDs[,,1]+hvD[,2]+hvDs[,,3]+hvD[,4]+hvDs[,,5]
		pvDs_rel <- vDs_rel / colSums(vDs_rel)
		E_Is_rel <- E_Is - round(colSums(pvDs_rel * I_NMDA * 1),3)
		# clustind <- round((-1) * E_Is / Rin * colSums(dpvDs * I_NMDA2),5)
		# main = paste('E[I]:', round(mean(E_Is), 2), 'dI/dn:', round(mean(clustind), 5))

		plot(vv, I_NMDA*ww, col=viridis(7)[5], lwd=2, t='l', xlab='', ylab='', xlim=c(-80, 10), ylim=c(-15, 20), axes=F, main=main); axis(2, las=2); abline(h=0)
		pvD <- rowMeans(pvDs)
		dpvD <- rowMeans(dpvDs)
		polygon(c(vv, rev(vv)), c(pvD*150, rep(0, length(vv))), col=grey(0.75, alpha=0.5))
		polygon(c(dvv, rev(dvv)), c(dpvD*300, rep(0, length(dvv))), col=grey(0.25, alpha=0.5))

		if (i.clust == 7){
			axis(1)
			mtext('voltage (mV)', 1, 2, cex=0.7)
			mtext('current (pA)', 2, 2.5, cex=0.7)
			# mtext('dI / dV (pA/mV)', 2, 2.5, cex=0.7)
		}		
				
		# clust_index[ww, i.clust,,] <- matrix(clustind, 10, 4)
		E_Vs <- colSums(pvDs * vv)
		E_Vd[ww, i.clust,,] <- matrix(E_Vs, 10, 4)
		E_I[ww, i.clust,,] <- E_Is
		E_I_rel[ww, i.clust,,] <- E_Is_rel
	}
}


################################
## plot the dendritic voltage
library(gplots)
setwd(localdir)

###################################################################################
### Fig4d, top, place / theta case
par(mar=c(4,4,1,1))
m1 <- apply(Vd[1,,], 1, mean); s1 <- apply(Vd[1,,], 1, sd) / sqrt(40)
m2 <- apply(Vd[2,,], 1, mean); s2 <- apply(Vd[2,,], 1, sd) / sqrt(40)
plotCI(1:7, m1, s1, pt.bg=viridis(8, option='D')[5], col=viridis(8, option='D')[3], t='o', pch=22, lty=1, ylim=c(-70, 0), axes=F, xlab='cluster size', ylab='Vm dendrite, mV', main='', gap=0)
plotCI(1:7, m2, s2, pt.bg=viridis(8, option='B')[5], col=viridis(8, option='B')[3], t='o', pch=22, lty=1, add=T, gap=0)
axis(2, las=2); axis(1, 1:7, nSynPerClust)

ampl_data <- rbind(Vd[1,,], Vd[2,,])
rnames <- c('240_cluster_of_1', '120_cluster_of_2', '48_cluster_of_5', '24_cluster_of_10', '12_cluster_of_20', '8_cluster_of_30', '4_cluster_of_60')
rownames(ampl_data) <- paste(rep(c('dVm_place', 'dVm_place_LTP'), e=7), rnames, sep='')
write.table(format(round(ampl_data, 3), scientific=F), file='Fig4d_topP_data.txt')

