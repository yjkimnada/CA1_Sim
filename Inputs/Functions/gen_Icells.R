## Functions to generate synthetic inhibitory cell activity in visual cortex
## simulates response to oriented moving gratings


#################################################################

gen.V1ICells <- function(N.cells, Tmax =24, Ntrial=16, dt=0.001, rate.ave=10, seed=23, graphics=T, adapt=T, gamma=NULL, oris=NULL){

# N.cells: number of cells to be simulated
# Tmax: length of the current maze, s
# Ntrial: number of trials simulated, 
# dt: temporal resolution (s)
# rate.ave: average firing rate
# seed: random seed
# graphics: if TRUE, results are plotted
# adapt: spike adaptation is modeled
# gamma: the phase of the gamma oscillation

	# N.cells <- Insyn; Ntrial <- 2; dt <- 0.001; seed = 52; rate.ave <- Irate; graphics <- T; adapt <- T; gamma=sp.E$gamma; oris <- sp.E$oris

	N.cells <- ceiling(N.cells / 16) * 16 # equal orientation...
	cat('We will generate spikes for ', N.cells, 'neuron. \n')
	
	if (is.null(gamma)){
		freq.gamma <- 40 # Hz
		gamma <- sim.run.circ(Tmax = Tmax, xmax=2*pi, dt=dt, mu.v=freq.gamma*2*pi, sd.v=100/1000, tau.v=0.5, seed=seed, v.max= freq.gamma*1.5*2*pi, v.min= freq.gamma*0.5*2*pi, init.v=NULL)$x
		L <- Tmax / dt
	} else {
		L <- length(gamma)
		Tmax <- L * dt
	}


	if (is.null(oris)){
		Noris <- 16
		L.ori <- L / Noris
		oris <- rep(1:Noris*2*pi / Noris, each=L.ori)
	} else {
		if (L != length(oris)) stop(paste('length of gamma and orientations should be the same, L_gamma=', L, ', L_oris=', length(oris), sep=''))
	}

	# generate time stamps for the trials
	tStartEnd <- matrix(NA, Ntrial, 2)
	for (i.trial in 1:Ntrial) {
		tStartEnd[i.trial,] <- c(0, Tmax) + (i.trial - 1) * Tmax
	}

	
	#################################################		
	## generate basis functions
	cat('generating basis functions... \n')
	N.ori.basis <- 16
	phi.ori <- seq(0,2*pi, length=N.ori.basis+1)
	# phi.ori <- seq(0,2*pi, length=10*N.ori.basis)
	ori.basis <- gen.basis.circ.mat(x=phi.ori, n=N.ori.basis, graphics=F)

	N.phase.basis <- 6
	phi.phase <- seq(0,2*pi, length=N.phase.basis+1)
	# phi.phase <- seq(0,2*pi, length=10*N.phase.basis)
	phase.basis <- gen.basis.circ.mat(x=phi.phase, n=N.phase.basis, graphics=F)

	all_phi.basis <- matrix(NA, N.ori.basis + N.phase.basis, N.ori.basis  * N.phase.basis)
	ii <- 1
	for (i.ori in 1: N.ori.basis){
		for (i.gamma in 1: N.phase.basis){
			all_phi.basis[,ii] <- c(ori.basis[,i.ori], phase.basis[,i.gamma])
			ii <- ii+ 1
		}
	}
	phi <- rbind(rep(1, N.ori.basis  * N.phase.basis), all_phi.basis)

	L.basis <- 100
	pphi <- seq(0,2*pi, length=L.basis)
	ori.basis <- gen.basis.circ.mat(x=pphi, n=N.ori.basis, F)
	phase.basis <- gen.basis.circ.mat(x=pphi, n=N.phase.basis, F)

	w.gamma <- c(-.5, 0.2, 0.5, 0.2, -.5, -2)
	# plot(exp(t(w.gamma %*% phase.basis)), t='l', lty=1)

	w.ori <- c(1,2,3,4,3,2,1,0,1,2,3,4,3,2,1,0)/10
	# plot(exp(t(w.ori %*% ori.basis)), t='l', lty=1)

	ratemap <- exp(as.vector(c(w.ori, w.gamma) %*% all_phi.basis))
	# plot(ratemap)
	mean.rate <- mean(ratemap)
	w0 <- log(rate.ave / mean.rate) 
	# ratemap <- exp(as.vector(w.gamma %*% phi.basis)) * exp(w0)
	# lines(ratemap, t='l', col=2)
	
	ww <- matrix(rep(c(w.ori, w.gamma), N.cells / N.ori.basis), N.cells/N.ori.basis, byrow=T)
	for (i in 2:N.ori.basis){
		ww2 <- matrix(rep(c(shift(w.ori, i), w.gamma), N.cells / N.ori.basis), N.cells/N.ori.basis, byrow=T)
		ww <- rbind(ww, ww2)
	}
	
	w0 <- rep(w0, N.cells)
	
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


	#############################################################
	## basis activations and firing rate in the function of time
	cat('basis activations... \n')

	phi.t <- matrix(NA, N.phase.basis+N.ori.basis, L)
	for (t in 1:L) {
		i.ori <- which.min((pphi - oris[t])^2)
		i.gam <- which.min((pphi - gamma[t])^2)
		phi.t[,t] <- c(ori.basis[,i.ori], phase.basis[,i.gam])
	}

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
	plotCells <- 1:min(N.cells, 240)
	NplotCells <- length(plotCells)
	
	cat('calculating spikes... \n')
	for (i.trial in 1:Ntrial){	
		# rate.t <- exp(w.x %*% x.t) * exp(w0)
		set.seed(seed*10000 + i.trial * 100)
		resp <- sim.glm(phi.t, ww, cons, c.basis, w.c, w.self, w0, dt)
		spt <- resp$sp
		
		if (graphics){
			raster <- matrix(0, length(plotCells), L)
			for (i in which(spt[,1] %in% plotCells)) {
				raster[spt[i,1], spt[i,2]] <- raster[spt[i,1], spt[i,2]] + 1
			}
			rasters[[i.trial]] <- raster
		}
		
		spiketimes[[i.trial]] <- spt
		spiketimes[[i.trial]][,2] <- spiketimes[[i.trial]][,2] * dt + tStartEnd[i.trial,1]
		cat ("trial", i.trial, "finished. \n")
	}
	
	
	if (graphics){

		# empirical orientation tuning of cells - oris
		# empirical phase tuning of cells - phases
		# empirical gamma tuning of cells - gamma


		allRasts <- rasters[[1]]; for (i in 2:Ntrial) allRasts <- cbind(allRasts, rasters[[i]])
		allGamma <- rep(gamma, Ntrial)
		allOris <- rep(oris, Ntrial)
		# sp <- raster2sp(allRrasts)
		# Tmax <- length(x) * dt
		
		pdf(file='V1_Iinput_tuning.pdf', 10, 5, useDingbats=F)
		layout(matrix(c(1,1,1,2,3,4), 2, byrow=T))
		par(mar=c(4,4,1,1))

		plot(spt[1001:2000,2]/1000, spt[1001:2000, 1], pch='|', xlab='time (s)', ylab='cells', axes=F, col=2, ylim=c(0, N.cells)); axis(1); axis(2, las=2)
		rr <- colSums(raster)
		rr.f <- filter(rr, rep(1/10, 10))
		lines(spt[1000,2]:spt[2000,2]/1000, rr.f[spt[1000,2]:spt[2000,2]])

		brOri <- seq(0, 2*pi, length=(Noris+1))+pi/Noris
		counts.ori <- hist(allOris, br=brOri, plot=F)$counts * dt
		mids.ori <- hist(allOris, br=brOri, plot=F)$mids
		rates.measured <- matrix(NA, NplotCells, Noris)
		for (i in 1:NplotCells) rates.measured[i,] <- hist(allOris[allRasts[i,] > 0], br=brOri, plot=F)$counts / counts.ori 
		# polar.plot(rates.measured, mids.ori/pi*180, rp.type='p', line.col=rainbow(NplotCells), main='orientation tuning')
		polar.plot(rates.measured[1:10,], mids.ori/pi*180, rp.type='p', line.col=rainbow(10), main='orientation tuning', radial.lim=c(0, max(rates.measured)), show.grid.labels=3)

		brPhases <- seq(0, 2*pi, length=18)
		counts.gamma <- hist(allGamma, br=brPhases, plot=F)$counts * dt
		mids.gamma <- hist(allGamma, br=brPhases, plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 17)
		for (i in 1:N.cells) rates.measured[i,] <- hist(allGamma[allRasts[i,] > 0], br=brPhases, plot=F)$counts / counts.gamma 
		polar.plot(rates.measured[1:10,], mids.gamma/pi*180, rp.type='p', line.col=viridis(10), main='gamma tuning', radial.lim=c(0, max(rates.measured)))
		dev.off()
	}
	
	data <- list(spt=spiketimes, gamma=gamma)
}
