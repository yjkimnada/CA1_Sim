## A script to analyse the effect of clustering on the postsynaptic Vm
##########################################

source('./CombLin.R', chdir = TRUE)
source('./Extract_spiketimes.R', chdir = TRUE)
source('./remove_spikes.R', chdir = TRUE)
source('./scalebar.R', chdir = TRUE)

library(viridis)
library(colormap)
library(gplots)


eval_Vm_stats <- function(n.sim, stim_types, act_types, typenames, graphics=0, gN='0.8', ramp=F, sampling_freq=1000){
	if (ramp) {
		sds <- array(NA, dim=c(6, n.sim, 10), dimnames=list(c("Sd_Vm_filt", "mean_ttt_variance_filt", "var_Vm_peak", "mean_ttt_variance_peak", "sd_ttt_variance_peak", "ramp area"), typenames, paste("rep", 1:10)))
	} else {
		sds <- array(NA, dim=c(5, n.sim, 10), dimnames=list(c("Sd_Vm_filt", "mean_ttt_variance_filt", "var_Vm_peak", "mean_ttt_variance_peak", "sd_ttt_variance_peak"), typenames, paste("rep", 1:10)))		
	}
	
	nSps <- array(0, dim=c(6, n.sim, 10), dimnames=list(c("m_pre", "m_in", "m_post", "sd_pre", "sd_in", "sd_post"), typenames, paste("rep", 1:10)))

	meanresps <- array(0, dim=c(n.sim, 10, 10000))
	ttt_variance <- array(0, dim=c(n.sim, 10, 10000))
	ttt_variance_peak <- array(NA, dim=c(79, n.sim, 10))
	ntrial <- 16
	Tmax <- 10
	
	# graphics <- F
	
	for (i_type in 1:n.sim){
		
		act_type <- act_types[i_type]
		stim_type <- stim_types[i_type]
		if (act_type == 'active_NMDAc') Ntype <- 'noDendNa_' else Ntype <- ''
		if (stim_type == 'random_NR') randW <- '_randW' else randW <- ''
		for (stimseed in 1:10){
			cat(stimseed, "started \n")
		
				infile <- paste('./', act_type, '/vdata_T10_Ne2000_gA0.6_tauA1_gN', gN, '_Ni200_gG0.7_gB1.2_', Ntype, 'Er1_Ir7.5_', stim_type, '_rep16_stimseed', stimseed, randW, '.bin', sep='')
				# infile <- paste('./', act_type, '/vdata_T10_Ne2000_gA1.2_tauA1_gN1.6_Ni200_gG1.7_gB2.9_', Ntype, 'Er1_Ir7.5_', stim_type, '_rep16_stimseed', stimseed, randW, '.bin', sep='')
			if (file.exists(infile)){
				con <- file(infile, "rb")
				dim <- readBin(con, "integer", 2)
				resp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
				close(con)
						
				if (max(resp) > -40){
					nSps.seed <- matrix(0, ntrial, 3)
					for (i in 1:ntrial){
						resp.i <- resp[i,]
						sps <- extract.spiketimes(resp.i, sampling_freq, graphics=0,  limite=-30)
						spt <- sps$st
						if (length(spt>0)){
							nSps.seed[i, 1] <- sum((spt < 3*sampling_freq) & (spt >= 0.5*sampling_freq))
							nSps.seed[i, 2] <- sum((spt < 5.750*sampling_freq) & (spt >= 3.250*sampling_freq))
							nSps.seed[i, 3] <- sum((spt < 9.500*sampling_freq) & (spt >= 7*sampling_freq))
							if (sampling_freq > 1000){
								sps$delay <- sps$delay + 0.6
								sps$dur <- sps$dur + 1.2
							}
							rr <- remove.spikes(resp.i, sps, dt=1000/sampling_freq)
							resp[i,] <- rr[1,]
						}								
					}
					nSps[1:3, i_type, stimseed] <- apply(nSps.seed, 2, mean)
					nSps[4:6, i_type, stimseed] <- apply(nSps.seed, 2, sd)									
				}
				#########################################################################################
				## somatic Vm - removing theta

				ii <- seq(0, (10*sampling_freq), length=10001)[-1]
				resp.1 <- resp[,ii]

				
				sdfilt <- 0.1 # s
				filt <- dnorm(seq(-4, 4, length=sdfilt * 2 * 1000 * 8 + 1))
				filt <- filt / sum(filt)
				resp.filt <- t(filter(t(resp.1), filt, circular = T))
	
				# tt <- seq(5/200, 5, by=5/200) * 1000
				# meas <- colMeans(matrix(colMeans(resp.filt), 25))
				
				if (graphics == TRUE){	
					matplot(t(resp.1), t="l", lty=1, col=grey(0.7,alpha=0.3), xlab="time (ms)", ylab="voltage (mV)", axes=F, ylim=c(-70, -60)); axis(1); axis(2, las=2)
					lines(colMeans(resp.1), col=1, lwd=2)
					
					matplot(t(resp.filt), t="l", lty=1, col=rgb(1, 0.6, 0.6), xlab="time (ms)", ylab="voltage (mV)", add=T)
					lines(colMeans(resp.filt), col=2, lwd=2)
		
					# lines(tt, linpred.raw, col=rgb(0.2, 1, 0.2), lwd=4, lty=1)
				
					# legend("topright", legend=c("output - raw", "output - filtered", "predicted"), col=c(1, 2, rgb(0.2, 1, 0.2)), lwd=c(3,3,4,2), lty=1)
					legend("topright", legend=c("output - raw", "output - filtered"), col=c(1, 2), lwd=c(3,3), lty=1)
				}
				meanresps[i_type, stimseed,] <- colMeans(resp.filt)
		
		
				# "Sd_Vm", "SE", "var_Vm_peak", "mean_ttt_variance_peak", "sd_ttt_variance_peak"
				x1 <- 500; x2 <- 9500
				sds[1, i_type, stimseed] <- sd(colMeans(resp.filt)[x1:x2])	 # was 500:5000
				
				sds[2, i_type, stimseed] <- mean(apply(resp.filt[,x1:x2], 2, var)) # mean(apply(resp.filt[,x1:x2], 2, sd)/4)		
				ttt_variance[i_type, stimseed, ] <- apply(resp.filt, 2, var)	
	
				ampls <- rep(NA, 125)
				# theta_cycles <- matrix(NA, 1264, 125)
				for (i.phase in 1:125){
					k <- seq(124, by=125, length=79)	+ i.phase 
					ampls[i.phase] <- mean(resp.1[,k])
					# theta_cycles[,i.phase] <- as.vector(resp.1[,k])	
					# matplot(t(resp.1[,k]), pch=1)
				}			
	
				i.phase <- which.max(ampls)
				k <- seq(124, by=125, length=79)	+ i.phase 
				if (graphics > 0) {
					matplot(t(resp.1[,k]), pch=2)
					lines(colMeans(resp.1[,k]), lwd=2)
				}
	
				sds[3, i_type, stimseed] <- var(colMeans(resp.1[,k]))
				sds[4, i_type, stimseed] <- mean(apply(resp.1[,k], 2, var) / 16)
				sds[5, i_type, stimseed] <- sd(apply(resp.1[,k], 2, var) / 16)
				ttt_variance_peak[, i_type, stimseed] <- apply(resp.1[,k], 2, var) / 16
				
				if (ramp){
					baseline <- rowMeans(resp.filt[,1:1000 + Tmax*1000 - 2000])	
					mean(rowSums(resp.filt - baseline) / 1000)
					sds[6, i_type, stimseed] <- mean(rowSums(resp.filt - baseline) / 1000)
				}				
				# hist(apply(resp.1[,k], 2, var) / 16)
				# abline(v=var(colMeans(resp.1[,k])), col=2)
			}
		}
	}
	
	stats_data <- list(sds=sds, meanresps=meanresps, ttt_variance=ttt_variance, ttt_variance_peak=ttt_variance_peak, nSps=nSps)
	stats_data
}

eval_replay_stats <- function(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=1000){

	ntrial <- 16
	nstart <- 20

	sds <- array(NA, dim=c(7, n.sim, 10), dimnames=list(c("SE_Vm_SPW", "SD_mean", "SD_max", "SD_fMean", "SD_fMax", 'control_Ampl', 'clust_Ampl'), typenames, paste("rep", 1:10)))
	nSps <- array(NA, dim=c(n.sim, 10, ntrial, nstart), dimnames=list(typenames, paste("seed", 1:10), paste("trial", 1:ntrial),  paste('tstart', 1:20)))
	meanresps <- array(NA, dim=c(n.sim, 10, 300, 20))
	
	for (i_type in 1:n.sim){
		act_type <- act_types[i_type]
		stim_type <- stim_types[i_type]
		cat("\n", act_type, "started \n")
			
		for (stimseed in 1:10){
			cat(stimseed, "started \n")
		
			infile <- paste('./', act_type, '/vdata_T0.3_Ne2240_gA0.6_tauA1_gN0.8_Ni200_gG0.7_gB1.2_Er9_Ir30_', stim_type, '_rep', ntrial, '_stimseed', stimseed, '.bin', sep="")
			con <- file(infile, "rb")
			dim <- readBin(con, "integer", 2)
			oresp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
			close(con)
			# matplot(t(oresp[,1000:2100]), t='l', lty=1, ylim=c(-75, -35), col=rainbow(16))
			#########################################################################################
			## somatic spikes
	
			resp <- oresp
			for (i in 1:ntrial){
				resp.i <- oresp[i,]
				Vm.mat <- matrix(resp.i, 1501)
				for (j in 1:nstart){
					Vm.ij <- Vm.mat[,j]
						
					sps <- extract.spiketimes(Vm.ij, 5000, graphics=0,  limite=-30, dvdt.th=5)
					spt <- sps$st
					# spt_list[[i]] <- spt
					nSps[i_type, stimseed, i, j] <- length(spt)
				}
				## removing spikes
				sps <- extract.spiketimes(oresp[i,], 5000, graphics=0,  limite=-30)
				if (length(sps$st>0)){
					sps$delay <- sps$delay + 1; 		sps$dur <- sps$dur + 1.2; 
					rr <- remove.spikes(oresp[i,], sps, dt=0.2)
					resp[i,] <- rr[1,]
				}
			}
			# matplot(t(resp[,1000:2100]), t='l', lty=1, ylim=c(-75, -35), col=1, add=T)
				
			#########################################################################################
			## somatic Vm - removing ripples
	
			# ii <- rep(121:200, 20) + rep(0:19, each=80) * 300
			# ij <- rep(1:120, 20) + rep(0:19, each=120) * 300
			# ik <- rep(201:300, 20) + rep(0:19, each=100) * 300
			# t <- seq(1, 6004)
			# matplot(t, t(resp), t='l', col=grey(0.7), lty=1)
			# resp2 <- resp; resp2[,ik] <- NA; resp2[,ij] <- NA
			# matplot(t, t(resp2), t='l', col=2, lty=1, add=T)
	
			ii.1kHz <- rep(seq(5, 1500, by=5), 20) + rep(0:19*1501, each=300)
			resp.1kHz <- resp[,ii.1kHz]
			sdfilt <- 0.0025 # s
			filt <- dnorm(seq(-4, 4, length=sdfilt * 2 * 1000 * 8 + 1))
			filt <- filt / sum(filt)
			resp.filt <- t(filter(t(resp.1kHz), filt, circular = T))
			
			if (graphics == TRUE){	
				# pdf(file=paste('resp_spines', stimseed, '_itrial', i, '_istart', j, '.pdf', sep=''), 4, 4)
				j <- 15; i <- 2
				kk <- ((j-1)*300+1):(j*300)
				kk.5kHz <- ((j-1)*1501+1):(j*1500)
				t.5kHz <- seq(0, by=0.2, length=30020)
				t.1kHz <- seq(1, by=1, length=6000)
	
				par(mar=c(2,2,2,2))
				matplot(t.5kHz[kk.5kHz], t(oresp[,kk.5kHz]), t='l', col=grey(0.7,alpha=1), lty=1, xlab='', ylab='', axes=F, ylim=c(-75, -30))
				matplot(t.1kHz[kk], t(resp.filt[,kk]), t='l', col=rgb(1, 0.6, 0.6), add=T, lty=1)
				lines(t.5kHz[kk.5kHz], oresp[i,kk.5kHz])
				lines(t.1kHz[kk], colMeans(resp.filt[,kk]), col=2, lwd=2, t='l')
				scalebar2(50, 5, '50 ms', '5 mV', 'topright')
				# dev.off()
		
			}
	
			# "SE_Vm_SPW", "SD_mean", "SD_max", "SD_fMean", "SD_fMax", 'mean_ampl', 'max_Ampl'
			ii.SPW <- rep(121:200, 20) + rep(0:19, each=80) * 300
			ii.baseline <- rep(21:100, 20) + rep(0:19, each=80) * 300
			sds[1, i_type, stimseed] <- mean(apply(resp.filt[,ii.SPW], 2, sd)) / sqrt(ntrial)		
	
			fVm.mat <- matrix(colMeans(resp.filt)[ii.SPW], 80)
			Vm.mat <- matrix(colMeans(resp.1kHz)[ii.SPW], 80)
			# matplot(Vm.mat, t='l', lty=1, col=rainbow(20, end=0.7))
			# matplot(fVm.mat, t='l', lty=1, col=rainbow(20, end=0.7), add=T)
			
			
			means.V.SPW <- apply(Vm.mat, 2, mean)
			sds[2, i_type, stimseed] <- sd(means.V.SPW)
	
			maxs.V.SPW <- apply(Vm.mat, 2, max)
			sds[3, i_type, stimseed] <- sd(maxs.V.SPW)
	
			means.fV.SPW <- apply(fVm.mat, 2, mean)
			sds[4, i_type, stimseed] <- sd(means.fV.SPW)
	
			maxs.fV.SPW <- apply(fVm.mat, 2, max)
			sds[5, i_type, stimseed] <- sd(maxs.fV.SPW)
	
			
			Vm.baseline.mat <- matrix(colMeans(resp.1kHz)[ii.baseline], 80)
			means.baseline <- apply(Vm.baseline.mat[61:80,], 2, mean)
			# matplot(Vm.baseline.mat, t='l', lty=1, col=rainbow(20, end=0.7))
			# abline(h=means.baseline, col=rainbow(20, end=0.7))
			
			sds[6, i_type, stimseed] <- mean((means.V.SPW - means.baseline)[1:5])
			sds[7, i_type, stimseed] <- mean((means.V.SPW - means.baseline)[10:11])
	
			meanresp <- matrix(colMeans(resp.filt), ncol=20)
			# meanresp <- matrix(colMeans(resp.1kHz), ncol=20)
			# matplot(meanresp, t='l', col=rainbow(30), lty=1)
			meanresps[i_type,stimseed,,] <- meanresp	
		}
	}

	stats_data <- list(sds=sds, meanresps=meanresps, nSps=nSps)
	stats_data		
}

eval_strong_stats <- function(n.sim, stim_types, act_types, typenames, graphics=0, sampling_freq=1000){
	### similar to eval_replay-stats, just adapted to the strong branches case, where we had only one synaptic configuration and 16 trials
	ntrial <- 16
	nstart <- 20

	ampl_clust <- array(NA, dim=c(n.sim, ntrial, nstart), dimnames=list(typenames, paste("trial", 1:ntrial),  paste('tstart', 1:20)))
	ampl_base <- array(NA, dim=c(n.sim, ntrial, nstart), dimnames=list(typenames, paste("trial", 1:ntrial),  paste('tstart', 1:20)))
	nSps <- array(NA, dim=c(n.sim, ntrial, nstart), dimnames=list(typenames, paste("trial", 1:ntrial),  paste('tstart', 1:20)))
	meanresps <- array(NA, dim=c(n.sim, 300, 20))
	
	for (i_type in 1:n.sim){
		act_type <- act_types[i_type]
		stim_type <- stim_types[i_type]
		cat("\n", act_type, "started \n")
					
		infile <- paste('./', act_type, '/vdata_T0.3_Ne2240_gA0.6_tauA1_gN0.8_Ni200_gG0.7_gB1.2_Er9_Ir30_', stim_type, '_rep', ntrial, '_stimseed1.bin', sep="")
		con <- file(infile, "rb")
		dim <- readBin(con, "integer", 2)
		oresp <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
		close(con)
		# matplot(t(oresp[,1000:2100]), t='l', lty=1, ylim=c(-75, -35), col=rainbow(16))
		#########################################################################################
		## somatic spikes

		resp <- oresp
		for (i in 1:ntrial){
			resp.i <- oresp[i,]
			Vm.mat <- matrix(resp.i, 1501)
			for (j in 1:nstart){
				Vm.ij <- Vm.mat[,j]
					
				sps <- extract.spiketimes(Vm.ij, 5000, graphics=0,  limite=-30, dvdt.th=5)
				spt <- sps$st
				# spt_list[[i]] <- spt
				nSps[i_type, i, j] <- length(spt)
			}
			## removing spikes
			sps <- extract.spiketimes(oresp[i,], 5000, graphics=0,  limite=-30)
			if (length(sps$st>0)){
				sps$delay <- sps$delay + 1; 		sps$dur <- sps$dur + 1.2; 
				rr <- remove.spikes(oresp[i,], sps, dt=0.2)
				resp[i,] <- rr[1,]
			}
		}
		# matplot(t(resp[,1000:2100]), t='l', lty=1, ylim=c(-75, -35), col=1, add=T)
			
		#########################################################################################
		## somatic Vm - removing ripples

		# ii <- rep(121:200, 20) + rep(0:19, each=80) * 300
		# ij <- rep(1:120, 20) + rep(0:19, each=120) * 300
		# ik <- rep(201:300, 20) + rep(0:19, each=100) * 300
		# t <- seq(1, 6004)
		# matplot(t, t(resp), t='l', col=grey(0.7), lty=1)
		# resp2 <- resp; resp2[,ik] <- NA; resp2[,ij] <- NA
		# matplot(t, t(resp2), t='l', col=2, lty=1, add=T)

		ii.1kHz <- rep(seq(5, 1500, by=5), 20) + rep(0:19*1501, each=300)
		resp.1kHz <- resp[,ii.1kHz]
		sdfilt <- 0.0025 # s
		filt <- dnorm(seq(-4, 4, length=sdfilt * 2 * 1000 * 8 + 1))
		filt <- filt / sum(filt)
		resp.filt <- t(filter(t(resp.1kHz), filt, circular = T))
		
		if (graphics == TRUE){	
			# pdf(file=paste('resp_spines', stimseed, '_itrial', i, '_istart', j, '.pdf', sep=''), 4, 4)
			j <- 15; i <- 2
			kk <- ((j-1)*300+1):(j*300)
			kk.5kHz <- ((j-1)*1501+1):(j*1500)
			t.5kHz <- seq(0, by=0.2, length=30020)
			t.1kHz <- seq(1, by=1, length=6000)

			par(mar=c(2,2,2,2))
			matplot(t.5kHz[kk.5kHz], t(oresp[,kk.5kHz]), t='l', col=grey(0.7,alpha=1), lty=1, xlab='', ylab='', axes=F, ylim=c(-75, -30))
			matplot(t.1kHz[kk], t(resp.filt[,kk]), t='l', col=rgb(1, 0.6, 0.6), add=T, lty=1)
			lines(t.5kHz[kk.5kHz], oresp[i,kk.5kHz])
			lines(t.1kHz[kk], colMeans(resp.filt[,kk]), col=2, lwd=2, t='l')
			scalebar2(50, 5, '50 ms', '5 mV', 'topright')
			# dev.off()
	
		}

		# "SE_Vm_SPW", "SD_mean", "SD_max", "SD_fMean", "SD_fMax", 'mean_ampl', 'max_Ampl'
		ii.SPW <- rep(121:200, 20) + rep(0:19, each=80) * 300
		ii.baseline <- rep(21:100, 20) + rep(0:19, each=80) * 300
		# sds[1, i_type, stimseed] <- mean(apply(resp.filt[,ii.SPW], 2, sd)) / sqrt(ntrial)		


		
		Vm.mat <- resp.1kHz[,ii.SPW]
		Vm.mat3 <- array(Vm.mat, dim=c(16, 80, 20))
		# plot(Vm.mat[3, 161:240])
		# lines(Vm.mat3[3, , 3])

		ampl_clust[i_type,,] <- apply(Vm.mat3, c(1,3), mean)

		Vm.base.mat <- resp.1kHz[,ii.baseline]
		Vm.base.mat3 <- array(Vm.base.mat, dim=c(16, 80, 20))
		ampl_base[i_type,,] <- apply(Vm.base.mat3[,61:80,], c(1,3), mean)


		meanresp <- matrix(colMeans(resp.filt), ncol=20)
		# meanresp <- matrix(colMeans(resp.1kHz), ncol=20)
		# matplot(meanresp, t='l', col=rainbow(30), lty=1)
		meanresps[i_type,,] <- meanresp	
	}

	stats_data <- list(ampl_clust = ampl_clust, ampl_base=ampl_base, meanresps=meanresps, nSps=nSps)
	stats_data		
}



plot_data <- function(filename=NULL, sds, ttt_variance, ttt_variance_peak, typenames, cols=NULL, pchs=NULL){
	if (!is.null(filename)) pdf(file=filename, 6, 3, useD=F)
	M <- dim(sds)[2]
	if (is.null(cols)) {
		icols <- seq(9, by=4, length=M)%%30 + 1
		cols <- colormap(colormaps$jet, nshades=30)[icols]
	}
	if (is.null(pchs)) {
		pchs <- rep(22, M)
	}
	par(mfcol=c(1,3))
	par(mar=c(6,4,2,1))
	# boxplot(t(sds[3,,]), ylim=c(0, 0.3), col=cols, axes=F, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)')
	
	mm <- apply(ttt_variance_peak, 2, mean)
	ss <- apply(ttt_variance_peak, 2, sd)
	ylim <- c(0, max(sds[3,,]))
	
	plotCI(1:M+0.15, mm, ss, gap=0, pch=17, lwd=1, ylim=ylim, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)', axes=F, xlim=c(0,M+1))
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	mm <- apply(sds[3,,], 1, mean)
	ss <- apply(sds[3,,], 1, sd)
	
	plotCI(1:M-0.15, mm, ss, gap=0, pch=pchs, add=T, lwd=1, pt.bg=cols, cex=2)
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	###################################
	mm <- apply(sds[1,,], 1, mean)
	ss <- apply(sds[1,,], 1, sd)
	ylim <- c(0, max(sds[1,,]))
	
	plotCI(1:M-0.15, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='SD of tuning', xlab='', ylab='voltage (mV)', axes=F, xlim=c(0,M+1), ylim=ylim)
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	# mm <- apply(sqrt(sds[2,,])/4, 1, mean)
	# ss <- apply(sqrt(sds[2,,])/4, 1, sd)
	mm <- apply(sqrt(ttt_variance)/4, 1, mean)
	ss <- apply(sqrt(ttt_variance)/4, 1, sd)
	
	plotCI(1:M+0.15, mm, ss, gap=0, pch=17, lwd=1, add=T)
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	#########################################
	mm <- apply(sds[1,,]^2, 1, mean)
	ss <- apply(sds[1,,]^2, 1, sd)
	ylim <- c(0, max(sds[1,,]^2))
	
	plotCI(1:M-0.15, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='variance of tuning', xlab='', ylab='voltage2 (mV2)', axes=F, xlim=c(0,M+1), ylim=ylim)
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	# mm <- apply(sds[2,,], 1, mean) / 16
	# ss <- apply(sds[2,,], 1, sd) / 16
	mm <- apply(ttt_variance, 1, mean) / 16
	ss <- apply(ttt_variance/16, 1, sd)
	
	plotCI(1:M+0.15, mm, ss, gap=0, pch=17, lwd=1, add=T)
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	if (!is.null(filename)) dev.off()
}




plot_data_SEM <- function(filename=NULL, sds, ttt_variance, ttt_variance_peak, typenames, cols=NULL, pchs=NULL, plotvars=F){
	if (!is.null(filename)) pdf(file=filename, 6, 3, useD=F)
	M <- dim(sds)[2]
	if (is.null(cols)) {
		icols <- seq(9, by=4, length=M)%%30 + 1
		cols <- colormap(colormaps$jet, nshades=30)[icols]
	}
	if (is.null(pchs)) {
		pchs <- rep(22, M)
	}
	par(mfcol=c(1,3))
	par(mar=c(6,4,2,1))
	# boxplot(t(sds[3,,]), ylim=c(0, 0.3), col=cols, axes=F, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)')
		
	mm <- apply(sds[3,,], 1, mean)
	ss <- apply(sds[3,,], 1, sd) / sqrt(10)
	ylim <- c(0, max(sds[3,,]+0.05))

	plotCI(1:M+0.15, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, xlim=c(0,M+1), ylim=ylim, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)', axes=F)
	matplot(1:M-0.15, sds[3,,], col=grey(0.74), pch=16, cex=0.7, add=T, t='p', lty=1)	
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	mm <- apply(ttt_variance_peak, 2, mean)	
	segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lwd=1, col=rgb(127, 42, 255, max=255))
	if (plotvars) {
		mm <- mm * 16	
		segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lty=3, lwd=1, col=rgb(255, 102, 0, max=255))
	}
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	
	###################################
	mm <- apply(sds[1,,], 1, mean)
	ss <- apply(sds[1,,], 1, sd) / sqrt(10)
	ylim <- c(0, max(sds[1,,]+0.05))
	
	plotCI(1:M+0.15, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='SD of tuning', xlab='', ylab='voltage (mV)', axes=F, xlim=c(0,M+1), ylim=ylim)
	matplot(1:M-0.15, sds[1,,], col=grey(0.74), pch=16, cex=0.7, add=T, t='p', lty=1)	
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	mm <- apply(sqrt(ttt_variance)/4, 1, mean)
	segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lwd=1, col=rgb(127, 42, 255, max=255))
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)

	if (plotvars) {
		mm <- mm * 4
		segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lty=3, lwd=1, col=rgb(255, 102, 0, max=255))
	}
	
	#########################################
	mm <- apply(sds[1,,]^2, 1, mean)
	ss <- apply(sds[1,,]^2, 1, sd) / sqrt(10)
	ylim <- c(0, max(sds[1,,]^2+0.02))
	
	plotCI(1:M+0.15, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='variance of tuning', xlab='', ylab='voltage2 (mV2)', axes=F, xlim=c(0,M+1), ylim=ylim)
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	matplot(1:M-0.15, sds[1,,]^2, col=grey(0.74), pch=16, cex=0.7, add=T, t='p', lty=1)	
	
	# mm <- apply(sds[2,,], 1, mean) / 16
	# ss <- apply(sds[2,,], 1, sd) / 16
	mm <- apply(ttt_variance, 1, mean) / 16
	segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lwd=1, col=rgb(127, 42, 255, max=255))
	axis(1, 1:M, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	if (plotvars) {
		mm <- mm * 16	
		segments(x0=1:M+0.25, x1=1:M-0.25, y0=mm, lty=3, lwd=1, col=rgb(255, 102, 0, max=255))
	}
	
	if (!is.null(filename)) dev.off()
}

plot_data_SEM_clust <- function(filename=NULL, sds, ttt_variance, ttt_variance_peak, typenames, cols=NULL, pchs=NULL, xx=NULL){
	if (!is.null(filename)) pdf(file=filename, 8, 2.5, useD=F)
	M <- dim(sds)[2]
	if (is.null(cols)) {
		icols <- seq(9, by=4, length=M)%%30 + 1
		cols <- colormap(colormaps$jet, nshades=30)[icols]
	}
	if (is.null(pchs)) {
		pchs <- rep(22, M)
	}
	if (is.null(xx)){
		xx <- 1:M
	}
	xlim = c(0, max(xx)+1)

	par(mfcol=c(1,4))
	par(mar=c(6,4,2,1))
	# boxplot(t(sds[3,,]), ylim=c(0, 0.3), col=cols, axes=F, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)')
		
	mm <- apply(sds[3,,], 1, mean)
	ss <- apply(sds[3,,], 1, sd) / sqrt(10)
	ylim <- c(0, 4)

	plotCI(xx, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, xlim=xlim, ylim=ylim, main='variability of theta peaks', xlab='', ylab='voltage2 (mV2)', axes=F)
	matplot(xx, sds[3,,], col=grey(0.74), pch=16, cex=0.7, add=T)	
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	mm <- apply(ttt_variance_peak, 2, mean)	
	segments(x0=xx+0.25, x1=xx-0.25, y0=mm, pch='-', lwd=1, axes=F, col=cols)
	axis(1, xx, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	
	###################################
	mm <- apply(sds[1,,], 1, mean)
	ss <- apply(sds[1,,], 1, sd) / sqrt(10)
	ylim <- c(0, 2)
	
	plotCI(xx, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='SD of tuning', xlab='', ylab='voltage (mV)', axes=F, xlim=xlim, ylim=ylim)
	matplot(xx, sds[1,,], col=grey(0.74), pch=16, cex=0.7, add=T)	
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	
	mm <- apply(sqrt(ttt_variance)/4, 1, mean)
	segments(x0=xx+0.25, x1=xx-0.25, y0=mm, pch='-', lwd=1, axes=F, col=cols)
	axis(1, xx, typenames, las=2, tick=F, line=-1); axis(2, las=2)
	
	#########################################
	mm <- apply(sds[1,,]^2, 1, mean)
	ss <- apply(sds[1,,]^2, 1, sd) / sqrt(10)
	ylim <- c(0, 3)
	
	plotCI(xx, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='variance of tuning', xlab='', ylab='voltage2 (mV2)', axes=F, xlim=xlim, ylim=ylim)
	legend('topleft', leg=c('measured', 'ttt expected'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
	matplot(xx, sds[1,,]^2, col=grey(0.74), pch=16, cex=0.7, add=T)	
	
	# mm <- apply(sds[2,,], 1, mean) / 16
	# ss <- apply(sds[2,,], 1, sd) / 16
	mm <- apply(ttt_variance, 1, mean) / 16
	segments(x0=xx+0.25, x1=xx-0.25, y0=mm, pch='-', lwd=1, axes=F, col=cols)
	axis(1, xx, typenames, las=2, tick=F, line=-1); axis(2, las=2)

	###################################
	mm <- apply(sds[6,,], 1, mean)
	ss <- apply(sds[6,,], 1, sd) / sqrt(10)
	ylim <- c(-2, 8.5)
	
	plotCI(xx, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, main='response integral', xlab='', ylab='mV x ms', axes=F, xlim=xlim, ylim=ylim)
	matplot(xx, sds[6,,], col=grey(0.74), pch=16, cex=0.7, add=T)	
	axis(1, xx, typenames, las=2, tick=F, line=-1); axis(2, las=2)

	if (!is.null(filename)) dev.off()
}


pplot <- function(data, baseline, xlab, cols=NULL, pchs=NULL, xx=NULL, ylim=NULL, add=F, greylevel=0.75, ...){
	M <- nrow(data)
	N <- ncol(data)
	if (is.null(cols)) {
		icols <- seq(9, by=4, length=M)%%30 + 1
		cols <- colormap(colormaps$jet, nshades=30)[icols]
	}
	if (is.null(pchs)) {
		pchs <- rep(22, M)
	}
	if (is.null(xx)){
		xx <- 1:M
	}
	if (is.null(ylim)){
		ylim = c(0, max(data))
	}

	#########################################
	mm <- apply(data, 1, mean)
	ss <- apply(data, 1, sd) / sqrt(N)
	
	plotCI(xx, mm, ss, gap=0, pch=pchs, lwd=1, pt.bg=cols, cex=2, xlab='', axes=F, ylim=ylim, add=add, ...)
	matplot(xx, data, bg=grey(greylevel), pch=21, cex=0.7, add=T, col=1)	
	if (add==F){	
		# legend('topleft', leg=c('measured', 'baseline'), pch=c(21, 24), pt.bg=c(cols[2], 1), bty='n')
		segments(x0=xx+0.25, x1=xx-0.25, y0=baseline, pch='-', lwd=1, axes=F, col=cols)
		axis(1, xx, xlab, las=1, tick=F, line=-1); axis(2, las=2)
	}
}
