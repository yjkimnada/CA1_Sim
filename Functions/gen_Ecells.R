## Functions to generate synthetic data that replicates visual cortical excitatory inputs activity
## simulates response to oriented moving gratings

source('./BasisFunctions.R', chdir = TRUE) # generating the basis functions
source('./SimRun.R', chdir = TRUE) # generating the basis functions
source("./sim_glm.R")
source('./raster2isp.R', chdir = TRUE)
# require('circular')
# require('plotrix')
# library('png')
# library('viridis')

##############################################################################
## location dependent cells - tuning to space, adaptation and velocity

shift <- function(x, n, circular=T, default=NA){
	m <- length(x)
	n <- n %% m

	if(n==0){
		return(x)
	}
	
	M <- length(x)
	if(n<0){
		n <- abs(n)
		forward=FALSE
	}else{
		forward=TRUE
	}
	if(forward){
		if (circular) return(c(x[seq(M-n+1, M)], x[seq(1, M-n)])) else return(c(rep(default, n), x[seq(1, M-n)]))
	} else {
		if (circular) return(c(x[seq(n+1, M)], x[seq(1, n)])) else  return(c(x[seq(n+1, M)], rep(default, n)))
	}
}



gen.V1Ecells <- function(Ns=c(ori=8, sdOri=3, phase=8, sdPhase=5, rep=2), Tmax =24, Ntrial=16, Noris=16, dt=0.001, rate.ave=1, seed=23, graphics=T, adapt=T, rand.tuning=T, rand.rates=F, oris=NULL, phases=NULL, gamma=NULL){

# Ns: number of cells to be simulated
# Tmax: length of the current maze, s
# Ntrial: number of trials simulated, 
# Noris: number of trials simulated, 
# dt: temporal resolution (s)
# rate.ave: average firing rate
# seed: random seed
# graphics: if TRUE, results are plotted
# adapt: spike adaptation is modeled
# rand.tuning: random tuning or uniform?
# rand.rates: random peak rate or the same

	# Ns <- c(ori=1, sdOri=3, phase=2, sdPhase=5); Tmax <- 20; Ntrial <- 10; Noris <- 8;  dt <- 0.001; seed = 52; N.basis <- 8; rate.ave <- 1; graphics <- T; adapt <- F; rand.tuning <- F; rand.rates <- F	
	# Ns <- ENs; seed = rseed; rate.ave <- Erate; graphics <- T; adapt <- F; rand.tuning <- rand.tuning; rand.rates <- rand.rates; phases <- NULL; oris <- NULL; gamma <- NULL	
	# Ns <- INs; Tmax <- Tmax; Ntrial <- Ntrial; Noris <- Noris;  dt <- 0.001; seed = rseed; rate.ave <- Irate; graphics <- T; adapt <- F; rand.tuning <- rand.tuning; rand.rates <- rand.rates; phases <- NULL; oris <- NULL; gamma <- NULL

	if (Ns[1] > 16) {
		Ns[1] <- 16
		warning('currently only 16 different orientations are defined')
	}
	if (Ns[2] > 5) {
		Ns[2] <- 5
		warning('currently only 5 different orientation tuning width are defined (orientation versus direction selectivity)')
	}
	if (Ns[3] > 6) {
		Ns[3] <- 6
		warning('currently only 6 different phase tuning can be defined')
	}
	if (Ns[4] > 4) {
		Ns[4] <- 4
		warning('currently only 4 different phase tuning width can be defined (simple and complex cells)')
	}

	N.cells <- prod(Ns)
	cat('We will generate spikes for ', N.cells, 'neuron. \n')
		
	## generate phases
	
	if (is.null(phases)){
		freq.phases <- 2 # Hz
		phases <- sim.run.circ(Tmax = Tmax, xmax=2*pi, dt=dt, mu.v=freq.phases*2*pi, sd.v=5/1000, tau.v=0.5, seed=seed, v.max=freq.phases*1.5*2*pi, v.min=freq.phases*0.5*2*pi, init.v=NULL)$x
		L <- Tmax / dt
	} else {
		cat('phases are provided \n')
		L <- length(phases)
		Tmax <- L * dt
	}

	if (is.null(gamma)){
		freq.gamma <- 40 # Hz
		gamma <- sim.run.circ(Tmax = Tmax, xmax=2*pi, dt=dt, mu.v=freq.gamma*2*pi, sd.v=100/1000, tau.v=0.5, seed=seed, v.max= freq.gamma*1.5*2*pi, v.min= freq.gamma*0.5*2*pi, init.v=NULL)$x
		L <- Tmax / dt
	} else {
		cat('gamma is provided \n')
		if (L != length(gamma)) stop(paste('length of gamma and phases should be the same, L_phases=', L, ', L_gamma=', length(gamma), sep=''))
	}

	if (is.null(oris)){
		Noris <- 16
		L.ori <- L / Noris
		oris <- rep(1:Noris*2*pi / Noris, each=L.ori)
		cat('generating ', Noris, 'discrete orientations \n')
	} else {
		cat('orientations are provided \n')
		if (L != length(oris)) stop(paste('length of gamma and orientations should be the same, L_gamma=', L, ', L_oris=', length(oris), sep=''))
	}
	
	# generate time stamps for the trials
	tStartEnd <- matrix(NA, Ntrial, 2)
	for (i.trial in 1:Ntrial) {
		tStartEnd[i.trial,] <- c(0, Tmax) + (i.trial - 1) * Tmax
	}

	#################################################		
	## generate basis functions
	N.ori.basis <- 16
	phi.ori <- seq(0,2*pi, length=N.ori.basis+1)
	# phi.ori <- seq(0,2*pi, length=10*N.ori.basis)
	ori.basis <- gen.basis.circ.mat(x=phi.ori, n=N.ori.basis, graphics=F)

	N.phase.basis <- 6
	phi.phase <- seq(0,2*pi, length=N.phase.basis+1)
	# phi.phase <- seq(0,2*pi, length=10*N.phase.basis)
	phase.basis <- gen.basis.circ.mat(x=phi.phase, n=N.phase.basis, graphics=F)

	# activation of the basis functions at different combinations
	all_phi.basis <- matrix(NA, N.ori.basis + N.phase.basis + N.phase.basis, N.ori.basis * N.phase.basis * N.phase.basis)
	ii <- 1
	for (i.ori in 1:N.ori.basis){
		for (i.phase in 1:N.phase.basis){
			for (i.gamma in 1: N.phase.basis){
				all_phi.basis[,ii] <- c(ori.basis[,i.ori], phase.basis[,i.phase], phase.basis[,i.gamma])
				ii <- ii+ 1
			}
		}
	}
	phi.allComb <- rbind(rep(1, N.ori.basis * N.phase.basis * N.phase.basis), all_phi.basis)

	L.basis <- 97
	pphi <- seq(0,2*pi, length=L.basis)
	ori.basis <- gen.basis.circ.mat(x=pphi, n=N.ori.basis, F)
	phase.basis <- gen.basis.circ.mat(x=pphi, n=N.phase.basis, F)
	
	## generate basis functions for connectivity
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=6))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }

	#################################################		
	## generate weigths for the  basis functions
	w.ori <- matrix(NA, N.cells, N.ori.basis)
	w.phase <- matrix(NA, N.cells, N.phase.basis)
	w0 <- rep(NA, N.cells)	

	ww.ori <- as.matrix(read.table('./Functions/ww.ori.txt'))
	# if (Ns['sdOri'] == 1) ww.ori[1,] <- c(1,2,3,4,3,2,1,0,1,2,3,4,3,2,1,0)/10
	# matplot(seq(0, 360, length=100), exp(t(ww.ori %*% ori.basis)), t='l', lty=1)
	# rowSums(exp(ww.ori %*% ori.basis))
	# abline(h=1)
	# abline(v=c(135, 175))

	ww.phase <- as.matrix(read.table('./Functions/ww.phase.txt'))
	ww.phase <- ww.phase[c(3,4,1,2),]
	# matplot(exp(t(ww.phase %*% phase.basis)), t='l', lty=1)
	# rowSums(exp(ww.phase %*% phase.basis))	

	w.gamma <- c(-.5, 0.2, 0.5, 0.2, -.5, -2)
	# plot(exp(as.vector(w.gamma %*% phase.basis)), t='l')
	# sum(exp(as.vector(w.gamma %*% phase.basis)))

	if (rand.tuning) set.seed(seed)
	
	ratemap <- matrix(NA, N.cells, N.ori.basis * N.phase.basis * N.phase.basis)
	
	for (i in 1:Ns['ori'])  {
		for (j in 1:Ns['phase'])  {				
			for (ii in 1:Ns['sdOri'])  {
				for (jj in 1:Ns['sdPhase'])  {
					
					# i <- 1; j <- 2; ii <- 2; jj <- 6
					i.cell <- jj + (i-1) * Ns['sdOri'] * Ns['phase'] * Ns['sdPhase'] + (j-1) * Ns['sdOri'] * Ns['sdPhase'] + (ii-1) * Ns['sdPhase']
					if (rand.tuning){
						ori <- ww.ori[ii,]
						shift.ori <- sample(1:N.ori.basis, 1)
						w.ori[i.cell,] <- shift(ori, shift.ori)
				
						ph <- ww.phase[jj,]
						shift.ph <- sample(1:N.phase.basis, 1)
						w.phase[i.cell,] <- shift(ph, shift.ph)
					} else {
						ori <- ww.ori[ii,]
						shift.ori <- i
						w.ori[i.cell,] <- shift(ori, shift.ori)
				
						ph <- ww.phase[jj,]
						shift.ph <- j
						w.phase[i.cell,] <- shift(ph, shift.ph)
						# cat(i.cell, ' ')
					}
					
					rmap <- exp(as.vector(c(w.ori[i.cell,], w.phase[i.cell,] , w.gamma) %*% all_phi.basis))
					mean.rate <- mean(rmap)
					if (rand.rates) r0 <- rgamma(1, 2*rate.ave, 2) else r0 <- rate.ave
					cat('ph:', shift.ph, 'ori:', shift.ori, 'r0', r0, '\n')
					w0[i.cell] <- log(r0 / mean.rate) 

					ratemap[i.cell,] <- exp(as.vector(c(w0[i.cell], w.ori[i.cell,], w.phase[i.cell,] , w.gamma) %*% phi.allComb))
					 # plot(ratemap[i.cell,], t='l')
				}
			}
		}
	}
	
	cat("w generated \n")

	# colMeans(matrix(colMeans(ratemap), 6))
	w.gam <- matrix(rep(w.gamma, N.cells), N.cells, byrow=T)

	# ww <- cbind(w0, w.ori, w.phase, w.gam)	
	# ratemap <- ww %*% phi.allComb
	# image(1:N.cells, 1:576, ratemap, col= viridis(24))

	# corvec <- cor(t(exp(ratemap)), exp(ratemap[1,] ))
	# plot(corvec)
	# abline(v=seq(1,8)*120, col=3)
	# cormat <- cor(t(exp(ratemap)))
	# image(cormat, col=viridis(24))

	#############################################################
	## basis activations and firing rate in the function of time
	ww <- cbind(w.ori, w.phase, w.gam)	

	phi.t <- matrix(NA, 2*N.phase.basis+N.ori.basis, L)
	for (t in 1:L) {
		i.ori <- which.min((pphi - oris[t])^2)
		i.phase <- which.min((pphi - phases[t])^2)
		i.gam <- which.min((pphi - gamma[t])^2)
		phi.t[,t] <- c(ori.basis[,i.ori], phase.basis[,i.phase], phase.basis[,i.gam])
	}
	# matplot(t(phi.t[1:16,]), t='l')
	
	# no connections
	N.wc <- 1
	w.c <- matrix(0, N.wc, 6)
	cons <- matrix(NA, 1, 2)
	
	# adaptation
	if (adapt){
		w.refr <- c(-4, -1, 0, 0, 0, 0)
		w.self <- w.refr
		w.self <- matrix(rep(w.self, N.cells), N.cells, 6, byrow=T)	
	} else {
		w.self <- matrix(0, N.cells, 6)
	}

	spiketimes <- list()
	rasters <- list()
	plotCells <- 1:min(N.cells, 24)
	NplotCells <- length(plotCells)
	
	for (i.trial in 1:Ntrial){	
		# rate.t <- exp(w.x %*% x.t) * exp(w0)
		set.seed(seed*10000 + i.trial * 100)
		resp <- sim.glm(phi.t, ww, cons, c.basis, w.c, w.self, w0, dt)
		spt <- resp$sp # in ms
		
		if (graphics){
			raster <- matrix(0, length(plotCells), L)
			for (i in which(spt[,1] %in% plotCells)) {
				raster[spt[i,1], spt[i,2]] <- raster[spt[i,1], spt[i,2]] + 1
			}
			rasters[[i.trial]] <- raster
		}
		
		spiketimes[[i.trial]] <- spt
		spiketimes[[i.trial]][,2] <- spiketimes[[i.trial]][,2] * dt + tStartEnd[i.trial,1] # in sec
		cat ("trial", i.trial, "finished. \n")
	}
	
	
	if (graphics){

		# empirical orientation tuning of cells - oris
		# empirical phase tuning of cells - phases
		# empirical gamma tuning of cells - gamma


		allRasts <- rasters[[1]]; for (i in 2:Ntrial) allRasts <- cbind(allRasts, rasters[[i]])
		allOris <- rep(oris, Ntrial)
		allPhases <- rep(phases, Ntrial)
		allGamma <- rep(gamma, Ntrial)
		# sp <- raster2sp(allRrasts)
		# Tmax <- length(x) * dt
		
		pdf(file='V1_Einput_tuning.pdf', 10, 5, useDingbats=F)
		layout(matrix(c(1,1,1,2,3,4), 2, byrow=T))
		par(mar=c(4,4,1,1))

		raster <- matrix(0, N.cells, L)
		for (i in 1:nrow(spt)) {
			raster[spt[i,1], spt[i,2]] <- raster[spt[i,1], spt[i,2]] + 1
		}
		rr <- colSums(raster)
		rr.f <- filter(rr, rep(1/10, 10))

		if (N.cells > 500){
			n1 <- 10000; n2 <- 50000
		} else {
			n1 <- 2500; n2 <- 5000			
		}
		plot(spt[n1:n2,2]/1000, spt[n1:n2, 1], pch='.', xlab='time (s)', ylab='cells', axes=F, col=2, ylim=c(0, N.cells)); axis(1); axis(2, las=2)
		lines(spt[n1,2]:spt[n2,2]/1000, 20*rr.f[spt[n1,2]:spt[n2,2]])


		brOri <- seq(0, 2*pi, length=(Noris+1))+pi/Noris
		counts.ori <- hist(allOris, br=brOri, plot=F)$counts * dt
		mids.ori <- hist(allOris, br=brOri, plot=F)$mids
		rates.measured <- matrix(NA, NplotCells, Noris)
		for (i in 1:NplotCells) rates.measured[i,] <- hist(allOris[allRasts[i,] > 0], br=brOri, plot=F)$counts / counts.ori 
		# polar.plot(rates.measured, mids.ori/pi*180, rp.type='p', line.col=rainbow(NplotCells), main='orientation tuning')
		if (N.cells > 500){
			polar.plot(rates.measured[1:20,], mids.ori/pi*180, rp.type='p', line.col=rep(rainbow(5), each=4), main='orientation tuning', radial.lim=c(0, max(rates.measured)), show.grid.labels=3)
		} else {
			polar.plot(rates.measured[1:NplotCells,], mids.ori/pi*180, rp.type='p', line.col=rep(rainbow(2), each=12), main='orientation tuning', radial.lim=c(0, max(rates.measured)), show.grid.labels=3) 
		}

		brPhases <- seq(0, 2*pi, length=19)
		counts.phases <- hist(allPhases, br=brPhases, plot=F)$counts * dt
		mids.phases <- hist(allPhases, br=brPhases, plot=F)$mids
		rates.measured <- matrix(NA, NplotCells, 18)
		for (i in 1:NplotCells) rates.measured[i,] <- hist(allPhases[allRasts[i,] > 0], br=brPhases, plot=F)$counts / counts.phases 
		if (N.cells > 500){
			polar.plot(rates.measured[1:20,], mids.phases/pi*180, rp.type='p', line.col=rep(rainbow(4), 5), main='phase tuning', radial.lim=c(0, max(rates.measured)))
		} else {
			polar.plot(rates.measured[1:NplotCells,], mids.phases/pi*180, rp.type='p', line.col=rep(rainbow(2), 12), main='phase tuning', radial.lim=c(0, max(rates.measured)))
		}


		brPhases <- seq(0, 2*pi, length=19)
		counts.gamma <- hist(allGamma, br=brPhases, plot=F)$counts * dt
		mids.gamma <- hist(allGamma, br=brPhases, plot=F)$mids
		rates.measured <- matrix(NA, NplotCells, 18)
		for (i in 1:NplotCells) rates.measured[i,] <- hist(allGamma[allRasts[i,] > 0], br=brPhases, plot=F)$counts / counts.gamma 
		polar.plot(rates.measured[1:NplotCells,], mids.gamma/pi*180, rp.type='p', line.col=viridis(24, option='B'), main='gamma tuning', radial.lim=c(0, max(rates.measured)))
		dev.off()
	}
	
	data <- list(spt=spiketimes, gamma=gamma, oris=oris, phases=phases, t=tStartEnd)
}

