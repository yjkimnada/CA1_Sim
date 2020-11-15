## Functions to generate synthetic data that replicates several aspacts of hippocampal activity

source('./BasisFunctions.R', chdir = TRUE) # generating the basis functions
source('./SimRun.R', chdir = TRUE) # generating the basis functions
source("./sim_glm.R")
source('./raster2isp.R', chdir = TRUE)

# library('png')
# library('viridis')

##############################################################################
## location dependent cells - tuning to space, adaptation and velocity
## can also generate sharp-wave like trajectories when sim.SPW=T with 15x faster speed and 5x higher firing rates
## and potentially larger arena

gen.locationcells <- function(N.cells=10, L =100, L.total=NULL, Ntrial=10, dt=0.01, N.x.basis=L.total/10+1, rate.ave=1, seed=23, mu.v=0.2, vx=NULL, phi=NULL, graphics=T, adapt=T, scale.speed=15, scale.rate=5, sim.SPW=F, rate.speed=F, multi.fields=F, sample.fields=F, randfields=F, randrates=F, fields.middle=F, var.trials=1, same.runs=F){

# N.cells: number of cells to be simulated
# L: length of the current maze
# L.total: length of the total arena encountered and replayed e.g., during SPWs
# Ntrial: number of trials simulated, 
# dt: temporal resolution (s)
# N.x.basis=L.total/10+1: number of spatial basis
# rate.ave: average firing rate
# seed: random seed
# vx, phi: possible position or theta phase vectors provided
# graphics: if TRUE, results are plotted
# adapt: spike adaptation is modeled
# scale.speed: factor of speed modulation during replay trajectories
# scale.rate: firing rate increase during replay
# sim.SPW: simulate sharp wave replay
# rate.speed: neuronal gain is modulated by the animal's speed
# multi.fields: random (T) fields, possibly multi-peaked as in Rich et al., 2014
# sample.fields: random (T) or stereoptye place fields
# fields.middle: place fields only at the middle of the track
# randfields: random location or uniform?
# randrates: random peak rate of the same
# same.runs: all runs are identical

# N.cells <- 10; L <- 100; L.total<- L; Ntrial <- 10; dt <- 0.001; seed = 52; N.x.basis <- 20; rate.ave <- 1; vx <- NULL; phi <- NULL; graphics <- T; adapt <- T; scale.speed=15; scale.rate=5, sim.SPW=T, rate.speed=F, multi.fields=F, sample.fields=F, randfields=T, randrates=T, fields.middle=F, var.trials=1, same.runs=F

# distance measured in cm!
	set.seed(seed)

	sd.v <- 0.1 * var.trials # standard deviation of OU-speed (default: 0.1)
	sd.phi <- 2.5 * var.trials # standard deviation of theta-speed (default: 2.5)
	
	
	print('trajectories started')


	if (is.null(L.total)) L.total <- L
	
	## generate a trajectory in 1D
	if (is.null(vx)){
		if (sim.SPW){
			vx <- sim.trial(N = Ntrial, xmax=L.total/100, dt =dt, seed=seed, sd.v = sd.v*scale.speed, mu.v=mu.v*scale.speed, tau.v=1, v.min=0.1*scale.speed, v.max=0.5*scale.speed, out.cm=T)
		} else {
			vx <- sim.trial(N = Ntrial, xmax=L/100, dt =dt, seed=seed, sd.v = sd.v, mu.v=mu.v, tau.v=1, v.min=0.1, v.max=0.5, out.cm=T)
			scale.rate <- 1
		}
		
		if (same.runs) {
			for (i in 2:length(vx$x)) {
				vx$x[[i]] <- vx$x[[1]]
				vx$v[[i]] <- vx$v[[2]]
			}
		}
	}
	
	
	print('SPW trajectories started')
	
	cat("movement and speed generated \n")
	
	if (L !=L.total) {
		if (L.total < L) cat("place fields do not cover the entire track! Length of maze:", L, "Length of place fields:", L.total, "\n")
		if (L.total > L) cat("some place fields are not explored during running! Length of maze:", L, "Length of place fields:", L.total, "\n")
	}
	
	# generate time stamps for the trials
	tStartEnd <- matrix(NA, Ntrial, 2)
	times <- list()
	t0 <- dt
	for (i.trial in 1:Ntrial) {
		ttimes <- seq(t0, by=dt, length=length(vx$x[[i.trial]]))
		times[[i.trial]] <- ttimes
		tStartEnd[i.trial,] <- range(ttimes)
		t0 <- ceiling(max(ttimes + 1))
	}

	## generate theta oscillation
	if (is.null(phi)){
		phi <- list()
		for (i.trial in 1:Ntrial){
			Tmax.i <- length(vx$x[[i.trial]]) * dt
			if (same.runs) sseed <- seed + 1 else sseed <- seed + 1 + i.trial
			pp <- sim.run.circ(Tmax = Tmax.i, xmax=2*pi, dt=dt, mu.v=8*2*pi, sd.v=sd.phi*2*pi, tau.v=0.5, seed=sseed, v.max=12*2*pi, v.min=5*2*pi, init.v=NULL)$x
			phi[[i.trial]] <- cbind(theta=NA*pp, amplitude=NA*pp, phase=pp)
		}
	}

	cat("theta generated \n")
	
	## generate basis functions for place cells
	xx <- seq(0,L.total)
	x.basis <- gen.Gauss.basis.mat(x=xx, sds=5, n=N.x.basis, circular=T) # was 10 before...
	# matplot(t(x.basis), t="l", col=rainbow(N.x.basis, end=0.7), lty=1, axes=F, xlab="distance (cm)", ylab="activation", main="spatial basis"); axis(1); axis(2, las=2)
	
	## generate basis functions for connectivity
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	N.c.basis <- 6
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }
		
	## generate weigths for the location -  basis functions
	w.x <- matrix(NA, N.cells, N.x.basis)
	w0 <- rep(NA, N.cells)
	w.template <- c(0,2, 2,1.5,1.8,2, 2.5,2.5, 0,0)
	nx <- length(w.template)
	
	if (fields.middle) max.shift <- N.x.basis - nx else max.shift <- N.x.basis
	if (randfields) i.shifts <- sample(seq(0, max.shift-1), N.cells, T) else i.shifts <- (1:N.cells) %% max.shift
	

	for (i.cell in 1:N.cells){
		if (multi.fields) {
			wLoc <- gen.wLoc(x.basis, L.total, wmax=c(4, 5), alpha=0.57, beta=1/0.14, rate.ave=rate.ave, active.cells = T, seed=317+i.cell)
			w.x[i.cell,] <- wLoc
			w0[i.cell] <- attr(wLoc, 'w0')
		} else {
			if (sample.template) w.temp <- runif(length(w.template), 0, 2*w.template) else w.temp <- w.template
			i.shift <- i.shifts[i.cell]
			if (randrates) r0 <- rgamma(1, 6*rate.ave, 6) else r0 <- rate.ave

			w.mat <- rep(0, N.x.basis+nx)
			w.mat[(1:nx)+i.shift] <- w.temp
			w.mat[1:nx] <- w.mat[1:nx] + w.mat[(1:nx)+N.x.basis]
			wLoc <- w.mat[1:N.x.basis]

			ratemap <- t(wLoc %*% x.basis)
			mean.rate <- mean(exp(ratemap))
			w0[i.cell] <- log(r0 / mean.rate)
			w.x[i.cell,] <- wLoc
		}
	}
	cat("w generated \n")
	# w.xx <- matrix(rexp(N.cells*N.x.basis, 1), N.cells, N.x.basis)
	# w.x <- sigm(w.xx, a=5, sl=5, th=2.3) # 90% around 0
	# # hist(w.x)
	
	# ## weight of the constant basis to scale the firing rate 
	# ratemap <- t(w.x %*% x.basis)
	# i.w.x <- sort(apply(ratemap, 2, which.max), ind=T)$ix
	# w.x <- w.x[i.w.x,]
	# ratemap <- t(w.x %*% x.basis)
	# matplot(ratemap, t="l", lty=1, col=rainbow(10))
	
	# mean.rate <- colMeans(exp(ratemap))
	# w0 <- log(rate.ave / mean.rate) 
	
	ratemap <- w.x %*% x.basis + w0
	i.w.x <- sort(apply(ratemap, 1, which.max), ind=T)$ix
	w.x <- w.x[i.w.x,]
	w0 <- w0[i.w.x]
	ratemap <- w.x %*% x.basis + w0
	plot(exp(ratemap[6,]), t="o")
	# matplot(t(exp(ratemap)), t="l", lty=1, col=rainbow(N.cells))
	ratemap.SPW <- w.x %*% x.basis + w0 + log(scale.rate)
	
	## basis activations and firing rate in the function of time
	raster <- list()
	spiketimes <- list()
	sp.cells <- rep(0, N.cells)
	tt.cells <- 0
	if (graphics) rate.t <- rep(0, N.cells)

	if (rate.speed) {
		w.x <- cbind(rep(1, N.cells), w.x)
	}
		
	for (i.trial in 1:Ntrial){
		TT <- length(vx$x[[i.trial]])
		x.t <- matrix(NA, N.x.basis, TT)
		for (t in 1:TT) {
			i.x <- which.min((xx - vx$x[[i.trial]][t])^2)
			x.t[,t] <- x.basis[,i.x]
		}
		# velocity tuning
		if (rate.speed) {
			# rate.x <- (vx$v[[i.trial]] - 25) / 72 # gives ~0.8, 1, 1.2 at v = [10, 25, 40] cm/s
			rate.x <- (vx$v[[i.trial]] - 25) / 20 # gives ~0.5, 1, 2 at v = [10, 25, 40] cm/s
			x.t <- rbind(rate.x, x.t)
		}
	
		# print(dim(w.x))
		# print(dim(x.t))
		if (graphics) rate.t <- cbind(rate.t, exp(w.x %*% x.t) * exp(w0))
		# rate.tt <- exp(w.x %*% x.t) * exp(w0)
		# matplot(t(rate.tt), t="l", lty=1, add=T, col=1)
					
		# no connections
		N.wc <- 1
		w.c <- matrix(0, N.wc, N.c.basis)
		cons <- matrix(NA, 1, 2)
		
		# adaptation
		if (adapt){
			if (N.c.basis != 6) stop("w.self should be provided")
			w.burst <- c(-4, 1, 2, 0, -1, -2)
			w.adapt <- c(-4, 0, 0, -1/8, -1/4, -1/4)
			w.regular <- c(-4, -2, -1, -1, -1, 0)		
			w.refr <- c(-4, -1, 0, 0, 0, 0)
			w.self <- w.refr
			w.self <- matrix(rep(w.self, N.cells), N.cells, N.c.basis, byrow=T)	
		} else {
			w.self <- matrix(0, N.cells, N.c.basis)
		}
		
		resp <- sim.glm(x.t, w.x, cons, c.basis, w.c, w.self, w0+log(scale.rate), dt)
		spt <- matrix(0, N.cells, TT)
		for (i in 1:nrow(resp$sp)) spt[resp$sp[i,1], resp$sp[i,2]] <- spt[resp$sp[i,1], resp$sp[i,2]] + 1
		raster[[i.trial]] <- spt
		spiketimes[[i.trial]] <- resp$sp[,c(2,1)]
		spiketimes[[i.trial]][,1] <- spiketimes[[i.trial]][,1] * dt + tStartEnd[i.trial,1]
		sp.cells <- sp.cells + rowSums(spt)
		tt.cells <- tt.cells + ncol(spt) * dt
		cat ("trial", i.trial, "finished. \n")
	}
	if (graphics) rate.t <- rate.t[,-1]

	acells <- matrix(1, N.cells, 4, dimnames=list(NULL, c("tetrode", "cell", "number of spikes", "firing rate")))
	acells[,2] <- 1:N.cells
	acells[,3] <- sp.cells
	acells[,4] <- sp.cells / tt.cells
	
	
	syn_run <- list(runtype=1, t=tStartEnd, times=times, dists=vx$x, speed.1D=vx$v, raster=raster, acells=acells, areas="CA1", spt=spiketimes) 
	if (!sim.SPW) syn_run$lfp <- list(phi)

	if (graphics){
		x <- unlist(vx$x)
		rasts <- raster[[1]]; for (i in 2:Ntrial) rasts <- cbind(rasts, raster[[i]])
		theta <- phi[[1]][,3]; for (i in 2:Ntrial) theta <- c(theta, phi[[i]][,3])
		sp <- raster2sp(rasts)
		Tmax <- length(x) * dt
		
		layout(matrix(c(1,4,2,3), 2), widths=c(1,4))
		par(mar=c(4,4,1,1))
		
		print("matplot")	
		counts.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$counts * dt
		mids.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(x[rasts[i,] > 0], br=seq(0, max(x), length=21), plot=F)$counts / counts.x 
		matplot(t(exp(ratemap[,1:(L+1)])), t="l", lty=1, col=rainbow(N.cells), lwd=2, axes=F, xlab="distance (cm)", ylab="rate (Hz)"); axis(1); axis(2, las=2)
		matplot(mids.x, t(rates.measured), t="l", lty=2, col=rainbow(N.cells), lwd=1, add=T)
				
		par(mar=c(1,3,1,4))
		maxsp <- min(1000, nrow(sp))
		plot(sp[1:maxsp,2], sp[1:maxsp,1], pch="|", axes=F, xlab="", ylab=""); axis(2, las=2); mtext("cells", 2, 2, cex=0.85)
		lines(x/10, t="l", col=2); mtext("position (cm)", 4, col=2, line=2); axis(4, c(0,50, 100)/12+1, c(0,50,100), col=2)
		par(mar=c(4,3,1,4))
		t <- seq(dt, Tmax, by=dt)
		
		maxt <- round(sp[maxsp,2])
		matplot(t[1:maxt], t(rate.t[,1:maxt]), t="l", lty=1, col=rainbow(N.cells), lwd=2, axes=F, xlab="time (s)", ylab=""); axis(2, las=2); axis(1); mtext("rate (Hz)", 2, 2, cex=0.85)

		counts.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$counts * dt
		mids.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(theta[rasts[i,] > 0], br=seq(0,2*pi, length=21), plot=F)$counts / counts.phi 
		par(mar=c(4,4,1,1))
		matplot(mids.phi, t(rates.measured), t="l", lty=1, col=rainbow(N.cells), lwd=1, axes=F, xlab="theta phase", ylab="rate (Hz)"); axis(1); axis(2, las=2)

	}
	
	syn_run
}

#################################################################################################
## theta-cells - tuning to theta + adaptation
gen.thetacells <- function(N.cells=10, L =100, Ntrial=10, dt=0.01, N.phi.basis=4, rate.ave=1, seed=23, vx=NULL, phi=NULL, graphics=T, adapt=T, same.cells=F, phi.template=NULL){
# N.cells: number of cells to be simulated
# L: length of the current maze
# Ntrial: number of trials simulated, 
# dt: temporal resolution (s)
# N.phi.basis: number of theta basis
# rate.ave: average firing rate
# seed: random seed
# vx, phi: possible position or theta phase vectors provided
# graphics: if TRUE, results are plotted
# adapt: spike adaptation is modeled
# same.cells: T or F, cells have identical theta-tuning

# N.cells <- 10; L <- 100; Ntrial <- 10; dt <- 0.01; seed = 52; N.phi.basis <- 4; rate.ave <- 2; adapt <- F
# N.cells <- 10; L <- 100; Ntrial <- 16; dt <- 0.001; seed = 52; N.phi.basis <- 4; rate.ave <- 16; adapt <- T; vx <- list(x=sp.E$dists, v=sp.E$speed.1D); phi <- sp.E$lfp[[1]]; graphics <- T


	set.seed(seed)
	## generate a trajectory in 1D
	if (is.null(vx)){
		vx <- sim.trial(N = Ntrial, xmax=L/100, dt =dt, seed=seed, sd.v = 0.1, mu.v=0.2, tau.v=1, v.min=0.1, v.max=0.5, out.cm=T)
	}

	# generate time stamps for the trials
	tStartEnd <- matrix(NA, Ntrial, 2)
	times <- list()
	t0 <- dt
	for (i.trial in 1:Ntrial) {
		ttimes <- seq(t0, by=dt, length=length(vx$x[[i.trial]]))
		times[[i.trial]] <- ttimes
		tStartEnd[i.trial,] <- range(ttimes)
		t0 <- ceiling(max(ttimes + 1))
	}

	## generate theta oscillation
	if (is.null(phi)){
		phi <- list()
		for (i.trial in 1:Ntrial){
			Tmax.i <- length(vx$x[[i.trial]]) * dt
			pp <- sim.run.circ(Tmax = Tmax.i, xmax=2*pi, dt=dt, mu.v=8*2*pi, sd.v=2.5*2*pi, tau.v=0.5, seed=seed+1, v.max=12*2*pi, v.min=5*2*pi, init.v=NULL)$x
			phi[[i.trial]] <- cbind(theta=NA*pp, amplitude=NA*pp, phase=pp)
		}
	}
	
	## generate basis functions for theta cells
	L.phi.basis <- 25
	pphi <- seq(0,2*pi, length=L.phi.basis)
	phi.basis <- gen.basis.circ.mat(x=pphi, n=N.phi.basis, F)
	# if (graphics) matplot(t(phi.basis), t="l", col=rainbow(12), lty=1, axes=F, xlab="distance (cm)", ylab="activation", main="spatial basis"); axis(1); axis(2, las=2)
	
	## generate basis functions for connectivity
	N.c.basis <- 6
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }
		
	## generate random weigths for the basis functions
	# ww.phi <- matrix(rexp(N.cells*N.phi.basis, 1), N.cells, N.phi.basis)
	# w.phi <- sigm(ww.phi, a=5, sl=5, th=2) # 75% around 0
	w.phi <- matrix(0, N.cells, N.phi.basis)
	if (is.null(phi.template)) phi.template <- 	c(0.9,1,1.2,1.4) #c(0,0.4,1.1,0.4)
	if (same.cells){
		# w.phi[,1:N.phi.basis] <- rep(c(1/2,1,1.6,2), each=N.cells)
		w.phi[,1:N.phi.basis] <- rep(phi.template, each=N.cells)
	} else {
		for (i.cell in 1:N.cells) {
			Nw <- rbinom(1, N.phi.basis, N.phi.basis / 6)
			if (Nw>0){
				iw <- sample(seq(1,N.phi.basis), 1)
				iw <- seq(iw, iw+Nw-1) %% N.phi.basis + 1
				w.cell <- rep(1, N.phi.basis)
				w.cell[iw] <- runif(Nw, 2, 5)
				w.phi[i.cell,] <- w.cell
			}
		}
	}
	
	# hist(w.phi)
	
	## weight of the constant basis to scale the firing rate 
	phimap <- t(w.phi %*% phi.basis)
	i.w.phi <- sort(apply(phimap, 2, which.max), ind=T)$ix
	w.phi <- w.phi[i.w.phi,]

	phimap <- t(w.phi %*% phi.basis)
	
	mean.rate <- colMeans(exp(phimap))
	if (same.cells) r0 <- rate.ave else r0 <- rgamma(N.cells, 2*rate.ave, 2)
	w0 <- log(r0 / mean.rate) 
	phimap <- t(t(phimap) + w0)	
	# matplot(exp(phimap), col=rainbow(N.cells), t="l", lty=1)
	# lines(rowMeans(exp(phimap)), lwd=3)
	
	## basis activations and firing rate in the function of time
	raster <- list()
	spiketimes <- list()
	sp.cells <- rep(0, N.cells)
	tt.cells <- 0
	if (graphics) rate.t <- rep(0, N.cells)

	for (i.trial in 1:Ntrial){
		TT <- length(vx$x[[i.trial]])
		phi.trial <- phi[[i.trial]][,3]
		phi.t <- matrix(NA, N.phi.basis, TT)
		for (t in 1:TT) {
			i.phi <- which.min((pphi - phi.trial[t])^2)
			phi.t[,t] <- phi.basis[,i.phi] 
		}
		if (graphics) rate.t <- cbind(rate.t, exp(w.phi %*% phi.t) * exp(w0))

		# no connections
		N.wc <- 1
		w.c <- matrix(0, N.wc, N.c.basis)
		cons <- matrix(NA, 1, 2)
		
		# adaptation
		if (adapt){
			if (N.c.basis != 6) stop("w.self should be provided")
			w.burst <- c(-4, 1, 2, 0, -1, 0)
			w.adapt <- c(-4, 0, 0, -1/8, -1/4, -1/4)
			w.regular <- c(-4, -2, -1, -1, -1, 0)		
			w.refr <- c(-4, -1, 0, 0, 0, 0)
			w.refr.inh <- c(-2, 0, 0, 0, 0, 0)
			w.self <- w.refr
			w.self <- matrix(rep(w.self, N.cells), N.cells, N.c.basis, byrow=T)	
		} else {
			w.self <- matrix(0, N.cells, N.c.basis)
		}

		resp <- sim.glm(phi.t, w.phi, cons, c.basis, w.c, w.self, w0, dt)
		spt <- matrix(0, N.cells, TT)
		for (i in 1:nrow(resp$sp)) spt[resp$sp[i,1], resp$sp[i,2]] <- spt[resp$sp[i,1], resp$sp[i,2]] + 1
		raster[[i.trial]] <- spt
		spiketimes[[i.trial]] <- resp$sp[,c(2,1)]
		spiketimes[[i.trial]][,1] <- spiketimes[[i.trial]][,1] * dt + tStartEnd[i.trial,1]
		sp.cells <- sp.cells + rowSums(spt)
		tt.cells <- tt.cells + ncol(spt) * dt
	}
	if (graphics) rate.t <- rate.t[,-1]

	acells <- matrix(1, N.cells, 4, dimnames=list(NULL, c("tetrode", "cell", "number of spikes", "firing rate")))
	acells[,2] <- 1:N.cells
	acells[,3] <- sp.cells
	acells[,4] <- sp.cells / tt.cells
	
	syn_run <- list(runtype=1, t=tStartEnd, times=times, dists=vx$x, speed.1D=vx$v, raster=raster, acells=acells, lfp=list(phi), areas="CA1", spt=spiketimes)

	if (graphics){
		x <- unlist(vx$x)
		rasts <- raster[[1]]; for (i in 2:Ntrial) rasts <- cbind(rasts, raster[[i]])
		theta <- phi[[1]][,3]; for (i in 2:Ntrial) theta <- c(theta, phi[[i]][,3])
		sp <- raster2sp(rasts)
		Tmax <- length(x) * dt
		
		layout(matrix(c(1,4,2,3), 2), widths=c(1,4))
		par(mar=c(4,4,1,1))
		
		print("matplot")
		counts.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$counts * dt
		mids.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(x[rasts[i,] > 0], br=seq(0, max(x), length=21), plot=F)$counts / counts.x 
		matplot(mids.x, t(rates.measured), t="l", lty=2, col=rainbow(N.cells), lwd=1, xlab="distance (cm)", ylab="rate (Hz)", axes=F); axis(1); axis(2, las=2)
				
		par(mar=c(1,3,1,4))
		t <- seq(dt, Tmax, by=dt); Tmm <- 200
		i.tt <- which(sp[,2] < Tmm)
		plot(sp[i.tt,2]*dt, sp[i.tt,1], pch="|", axes=F, xlab="", ylab="", ylim=c(0, N.cells), xlim=c(0, Tmm/1000)); axis(2, las=2); mtext("cells", 2, 2, cex=0.85)
		lines(t[1:Tmm], x[1:Tmm]/10, t="l", col=2); mtext("position (cm)", 4, col=2, line=2); axis(4, c(0,50, 100)/12+1, c(0,50,100), col=2); axis(1)
		par(mar=c(4,3,1,4))
		matplot(t[1:Tmm], t(rate.t[,1:Tmm]), t="l", lty=1, col=rainbow(N.cells), lwd=2, axes=F, xlab="time (s)", ylab=""); axis(2, las=2); axis(1); mtext("rate (Hz)", 2, 2, cex=0.85)
		
		
		par(mar=c(4,4,1,1))
		matplot(seq(0, 2*pi, length=L.phi.basis), exp(phimap), t="l", lty=1, col=rainbow(N.cells), lwd=2, axes=F, xlab="theta phase", ylab="rate (Hz)"); axis(1); axis(2, las=2)
		counts.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$counts * dt
		mids.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(theta[rasts[i,] > 0], br=seq(0,2*pi, length=21), plot=F)$counts / counts.phi 
		matplot(mids.phi, t(rates.measured), t="l", lty=2, col=rainbow(N.cells), lwd=1, add=T)
	}	
	syn_run
}


########################################################################
## place cells - that have place field, are theta-modulated and show phase precession
## can also simulate adaptation

gen.placecells <- function(N.cells=10, L =100, L.total=L, Ntrial=10, dt=0.01, N.x.basis=L.total/5, N.phi.basis=8, rate.ave=1, rate.bg=1/100, seed=23, mu.v=0.2, vx=NULL, phi=NULL, graphics=T, plot.phase.precess=T, adapt=T, sample.template=F, randfields=T, randrates=T, fields.middle=T, var.trials=1, same.runs=F, w.template=NULL){
# N.cells: number of cells to be simulated
# L: length of the current maze (cm)
# L.total: length of the total arena encountered and replayed e.g., during SPWs
# Ntrial: number of trials simulated, 
# dt: temporal resolution (s)
# N.x.basis=L.total/10+1: number of spatial basis
# N.phi.basis: number of theta basis
# rate.ave: average firing rate
# seed: random seed
# vx, phi: possible position or theta phase vectors provided
# mu.v: mean velocity in m/s (default: 20 cm/s)
# graphics: if TRUE, results are plotted
# plot.phase.precess: 
# adapt: spike adaptation is modeled
# sample.template # to generate phase precession either sample randomly from a template location-phase tuning matrix or use the template
# randfields: place fields are either distributed uniformly, on the entire track, or randomly, 
# fields.middle: T or F fields should be located only in the middle of the track
# var.trials: real, variability of the speed and theta within a given trial
# same.runs: T or F whether speed and theta shoulf be the same on different trials

# N.cells <- 2000; L  <- 200; L.total <- L; Ntrial <- 16; dt <- 0.001; N.x.basis <- 40; N.phi.basis <- 4; rate.ave <- 1/2; rate.bg = 1/100; seed <- 1; vx <- NULL; phi <- NULL; graphics <- T; adapt <- T; sample.template=F; randfields=T; randrates=T; fields.middle=F; var.trials=1/1000; same.runs=T

# N.cells <- 10; L  <- 200; L.total <- L; Ntrial <- 10; dt <- 0.001; N.x.basis <- 40; N.phi.basis <- 4; rate.ave <- 1/2; rate.bg = 1; seed <- 1; vx <- NULL; phi <- NULL; graphics <- T; adapt <- T; sample.template=F; randfields=T; randrates=T; fields.middle=F; var.trials=1/100000; same.runs=T; w.template=w.template; mu.v=0.2

	set.seed(seed)
	if (is.null(w.template)) stop("w.template must be provided!")	
	if (ncol(w.template) != N.phi.basis) stop("ncol(w.template) must equal to N.phi.basis")
	nx <- nrow(w.template)
	
	# if (N.phi.basis !=5) stop("number of phi basis should be 5")
	# if (N.x.basis < 20) stop("number of x basis should be larger than 20")

	cat("generating phase precession... \n")
	cat("template sampled: ", sample.template, "; field location random: ", randfields, "; fields on the middle: ", fields.middle, "\n", sep="")

	sd.v <- 0.1 * var.trials # standard deviation of OU-speed (default: 0.1)
	sd.phi <- 2.5 * var.trials # standard deviation of theta-speed (default: 2.5)

	cat("sd of speed: ", sd.v, "; sd.phi: ", sd.phi, "\n", sep="")
	
	## generate a trajectory in 1D
	if (is.null(vx)){
		vx <- sim.trial(N = Ntrial, xmax=L/100, dt =dt, seed=seed, sd.v = sd.v, mu.v=mu.v, tau.v=1, v.min=0.1, v.max=2.5*mu.v, out.cm=T)
		if (same.runs) {
			for (i in 2:length(vx$x)) {
				vx$x[[i]] <- vx$x[[1]]
				vx$v[[i]] <- vx$v[[2]]
			}
		}
	}
		
	tStartEnd <- matrix(NA, Ntrial, 2)
	times <- list()
	t0 <- dt
	for (i.trial in 1:Ntrial) {
		ttimes <- seq(t0, by=dt, length=length(vx$x[[i.trial]]))
		times[[i.trial]] <- ttimes
		tStartEnd[i.trial,] <- range(ttimes)
		t0 <- ceiling(max(ttimes + 1))
	}
	
	## generate theta oscillation
	if (is.null(phi)){
		phi <- list()
		for (i.trial in 1:Ntrial){
			Tmax.i <- length(vx$x[[i.trial]]) * dt
			if (same.runs) sseed <- seed + 1 else sseed <- seed + 1 + i.trial
			pp <- sim.run.circ(Tmax = Tmax.i, xmax=2*pi, dt=dt, mu.v=8*2*pi, sd.v=sd.phi*2*pi, tau.v=0.5, seed=sseed, v.max=12*2*pi, v.min=5*2*pi, init.v=NULL)$x
			phi[[i.trial]] <- cbind(theta=NA*pp, amplitude=NA*pp, phase=pp)
		}
	}
	
	##################################################
	## 2. generate basis functions
	N.basis <- N.x.basis * N.phi.basis
	
	## 2D basis functions:
	## position, theta, direction
	ranges <- list(c(0, L.total), c(0, 2*pi)); ns <- c(N.x.basis, N.phi.basis)
	basis <- gen.Gauss.basis(ranges, ns, circular=c(T,T), force.sd=c(L.total/N.x.basis, NA))
	
	## generate basis functions for connectivity
	N.c.basis <- 6
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }
	
	# image(w.template, col=viridis(25))	

	set.seed(seed)
	ww <- matrix(NA, N.cells, N.x.basis * N.phi.basis)
	if (fields.middle) max.shift <- N.x.basis - nx else max.shift <- N.x.basis
	if (randfields) i.shifts <- sample(seq(0, max.shift-1), N.cells, T) else i.shifts <- (1:N.cells) %% max.shift
	
	for (i.cell in 1:N.cells){
		if (sample.template) w.temp <- matrix(runif(length(w.template), 0, 2*w.template), ncol=5) else w.temp <- w.template
		i.shift <- i.shifts[i.cell]
		w.mat <- matrix(0, N.x.basis + nx, N.phi.basis)
		w.mat[(1:nx)+i.shift,] <- w.temp
		w.mat[1:nx,] <- w.mat[1:nx,] + w.mat[(1:nx)+N.x.basis,]
		w.mat <- w.mat[1:N.x.basis,]
		# print(i.shift)
		# image(w.mat, col=viridis(25))	
		ww[i.cell,] <- as.vector(w.mat)
	}	
	w.x <- ww
	# unique(rowSums(w.x)) # should be length 1 if sample.template == F
	
	## weight of the constant basis to scale the firing rate 
	## the average firing rate of all cells is rate.ave (2 Hz)
	## r = f(x) g(theta) exp(w0)
	## E[r] = E[f] E[g] exp(w0)
	# source('../Code/glm4hpc/functs_hpc/RateXTS.R')
	
	Nx <- 100; Nphi <- 10
	xx <- seq(0,L.total, length=Nx+1)
	pphi <- seq(0,2*pi, length=Nphi+1)
	
	rates <- array(NA, dim=c(N.cells, Nx+1, Nphi+1))
	xphi.basis <- array(NA, dim=c(Nx+1, Nphi+1, N.basis))
	for (i.x in 1:(Nx+1)){
		x <- xx[i.x]
		for (i.phi in 1:(Nphi+1)){
			p <- pphi[i.phi]
			xphi.basis.i <- eval.Gauss.basis(c(x, p), basis, circular=T, period=c(L, 2*pi))
			rates[,i.x, i.phi] <- exp(w.x %*% xphi.basis.i)
			xphi.basis[i.x,i.phi, ] <- xphi.basis.i
		}
	}

	# for (i in 1:10) image(xx, pphi, rates[i,,], col=viridis(25))

	mean.rates <- apply(rates[,-1,-1], 1, mean)
	if (randrates) r0 <- rgamma(N.cells, 6*rate.ave, 6) else r0 <- rate.ave
	w0 <- log(r0 / mean.rates)
	
	ratemap <- apply(rates[,-1,-1], c(1,2), mean) * exp(w0)
	i.cell <- sort(apply(ratemap, 1, which.max), ind=T)$ix
	ratemap <- ratemap[i.cell,]
	w.x <- w.x[i.cell,]
	w0 <- w0[i.cell]

	# matplot(t(ratemap), t="l", col=rainbow(N.cells), lty=1)
	# plot(ratemap[4,], t="o")
	## basis activations and firing rate in the function of time
	raster <- list()
	spiketimes <- list()
	sp.cells <- rep(0, N.cells)
	tt.cells <- 0
	if (graphics) rate.t <- rep(0, N.cells)
	
	for (i.trial in 1:Ntrial){
		
		cat('starting trial',  i.trial, ' \n')
		cat('phi, ')
		phi.trial <- phi[[i.trial]][,3]
		cat('TT, ')
		TT <- length(vx$x[[i.trial]])
		cat('x.t with TT:', TT, '\n')
		x.t <- matrix(NA, N.basis, TT)

		for (t in 1:TT) {
			i.x <- which.min((xx - vx$x[[i.trial]][t])^2)
			i.phi <- which.min((pphi - phi.trial[t])^2)
			x.t[,t] <- xphi.basis[i.x,i.phi,]
		}
		if (graphics) rate.t <- cbind(rate.t, exp(w.x %*% x.t) * exp(w0))
		cat('basis functions calculated \n')
		
		# no connections
		N.wc <- 1
		w.c <- matrix(0, N.wc, N.c.basis)
		cons <- matrix(NA, 1, 2)
		
		# adaptation
		if (adapt){
			if (N.c.basis != 6) stop("w.self should be provided")
			w.burst <- c(-4, 1, 2, 0, -1, 0)
			w.adapt <- c(-4, 0, 0, -1/8, -1/4, -1/4)
			w.regular <- c(-4, -2, -1, -1, -1, 0)		
			w.refr <- c(-4, -1, 0, 0, 0, 0)
			w.self <- w.refr
			w.self <- matrix(rep(w.self, N.cells), N.cells, N.c.basis, byrow=T)	
		} else {
			w.self <- matrix(0, N.cells, N.c.basis)
		}
		
		cat('starting the simulations for trial',  i.trial, ' \n')
		set.seed(seed * 1000 + i.trial)
		resp <- sim.glm(x.t, w.x, cons, c.basis, w.c, w.self, w0, dt, rate.bg=rate.bg)
		cat("collecting the spikes \n")
		spt <- Matrix(0, N.cells, TT)
		for (i in 1:nrow(resp$sp)) spt[resp$sp[i,1], resp$sp[i,2]] <- spt[resp$sp[i,1], resp$sp[i,2]] + 1
		cat('spike raster prepared... \n')
		raster[[i.trial]] <- spt
		spiketimes[[i.trial]] <- resp$sp[,c(2,1)]
		spiketimes[[i.trial]][,1] <- spiketimes[[i.trial]][,1] * dt + tStartEnd[i.trial,1]
		sp.cells <- sp.cells + rowSums(spt)
		tt.cells <- tt.cells + ncol(spt) * dt
		cat("\n trial", i.trial, " finished. \n")
	}
	
	print("simulations finished")
	if (graphics) rate.t <- rate.t[,-1]

	acells <- matrix(1, N.cells, 4, dimnames=list(NULL, c("tetrode", "cell", "number of spikes", "firing rate")))
	acells[,2] <- 1:N.cells
	acells[,3] <- sp.cells
	acells[,4] <- sp.cells / tt.cells
	
	
	syn_run <- list(runtype=1, t=tStartEnd, times=times, dists=vx$x, speed.1D=vx$v, raster=raster, acells=acells, lfp=list(phi), areas="CA1", spt=spiketimes)

	############################

	cat("graphics:", graphics, "\n")
	
	if (graphics){
		x <- unlist(vx$x)
		rasts <- raster[[1]]; for (i in 2:Ntrial) rasts <- cbind(rasts, raster[[i]])
		theta <- phi[[1]][,3]; for (i in 2:Ntrial) theta <- c(theta, phi[[i]][,3])
		sp <- raster2sp(rasts)
		Tmax <- length(x) * dt
		
		layout(matrix(c(1,4,2,3), 2), widths=c(1,4))
		par(mar=c(4,4,1,1))
		
		counts.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$counts * dt
		mids.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(x[rasts[i,] > 0], br=seq(0, max(x), length=21), plot=F)$counts / counts.x 
		matplot(mids.x, t(rates.measured), t="l", lty=2, col=rainbow(N.cells, end=0.7), lwd=1, xlab="distance (cm)", ylab="rate (Hz)", axes=F); axis(1); axis(2, las=2)
		matplot(xx[-1], t(ratemap), t="l", col=rainbow(N.cells, end=0.7), lty=1, add=T)
		lines(xx[-1], colMeans(ratemap), t="l", lwd=2)

		par(mar=c(1,3,1,4))
		plot(sp[,2], sp[,1], xlim=c(0, 10000), pch="|", axes=F, xlab="", ylab=""); axis(2, las=2); mtext("cells", 2, 2, cex=0.85)
		lines(x*N.cells/L.total, t="l", col=2); mtext("position (cm)", 4, col=2, line=2); axis(4, c(0,50, 100)/12+1, c(0,50,100), col=2)

		par(mar=c(4,3,1,4))
		t <- seq(dt, Tmax, by=dt)
		matplot(t[1:10000], t(rate.t[, 1:10000] + rate.bg), t="l", lty=1, col=rainbow(N.cells, end=0.7), lwd=1, axes=F, xlab="time (s)", ylab=""); axis(2, las=2); axis(1); mtext("rate (Hz)", 2, 2, cex=0.85)
		
		counts.phi <- hist(theta, br=seq(0,2*pi, length=11), plot=F)$counts * dt
		mids.phi <- hist(theta, br=seq(0,2*pi, length=11), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 10)
		for (i in 1:N.cells) rates.measured[i,] <- hist(theta[rasts[i,] > 0], br=seq(0,2*pi, length=11), plot=F)$counts / counts.phi 
		par(mar=c(4,4,1,1))
		matplot(mids.phi, t(rates.measured), t="l", lty=2, col=rainbow(N.cells, end=0.7), lwd=1, axes=F, xlab="theta phase", ylab="rate (Hz)"); axis(1); axis(2, las=2)
			
		if (plot.phase.precess){
			sp <- cbind(sp, x[sp[,2]], theta[sp[,2]])
			sp[sp[,3] > L.total,3] <- L.total
			theta.time <- hist(theta, br=seq(0,2*pi, length=11), plot=F)$counts*dt
			cols <- c(rgb(1, 153/255, 85/255), rgb(42/255,127/255, 1))
			par(mfcol=c(5,2)); par(mar=c(3,3,1,1))
			for (i.cell in 1:10){
				isp.cell <- which(sp[,1]==i.cell)
				
				plot(sp[isp.cell, 3], sp[isp.cell, 4], pch=16, cex=0.7, axes=F, xlab="", ylab="", xlim=c(0,L), ylim=c(0,4*pi), col=rainbow(12)[i.cell])
				points(sp[isp.cell, 3], sp[isp.cell, 4]+2*pi, pch=16, cex=0.7, col=rainbow(12)[i.cell])
				placefield <- hist(sp[isp.cell, 3], br=seq(0, L.total, length=26), plot=F)
				lines(placefield$mids, placefield$counts/2/rate.ave, lwd=2, col=1)
						
				thetafield <- hist(sp[isp.cell, 4], br=seq(0,2*pi, length=11), plot=F)
				lines(c(thetafield$counts/theta.time, thetafield$counts/theta.time)*4, c(thetafield$mids, thetafield$mids+2*pi), lwd=2, col=rainbow(12)[i.cell])
			}
			axis(1, c(0, 50, 100), c(0, 50, 100))
			axis(2, c(0, 2*pi, 4*pi), c(0, expression(2* pi), expression(4 *pi)), las=2)
			mtext("postition (cm)", 1, 2, cex=0.7)
			mtext("theta phase", 2, 2, cex=0.7)
		}
	}
	syn_run

}



########################################################################
## place cells - that have place field, are ripple-modulated 
## can also simulate adaptation
## clustered cells do not show background activity - the background of their peers is separated into another number of neurons, put at the end of the block

gen.replay <- function(N.cells=10, L =100, L.total=L, Ntrial=10, dt=0.01, N.x.basis=L.total/5, N.phi.basis=8, rate.ave=50, rate.bg=0.4, seed=23, mu.v=5, vx=NULL, phi=NULL, graphics=T, plot.phase.precess=T, adapt=T, sample.template=F, randfields=T, randrates=T, fields.middle=T, var.trials=1, same.runs=F, w.template=NULL, phi.template=NULL, i.clust=NULL){
# N.cells: number of cells to be simulated
# L: length of the current maze (cm)
# L.total: length of the total arena encountered and replayed e.g., during SPWs
# Ntrial: number of trials simulated, 
# dt: temporal resolution (s)
# N.x.basis=L.total/10+1: number of spatial basis
# N.phi.basis: number of theta basis
# rate.ave: average firing rate
# seed: random seed
# vx, phi: possible position or theta phase vectors provided
# mu.v: mean velocity in m/s (default: 20 cm/s)
# graphics: if TRUE, results are plotted
# plot.phase.precess: 
# adapt: spike adaptation is modeled
# sample.template # to generate phase precession either sample randomly from a template location-phase tuning matrix or use the template
# randfields: place fields are either distributed uniformly, on the entire track, or randomly, 
# fields.middle: T or F fields should be located only in the middle of the track
# var.trials: real, variability of the speed and theta within a given trial
# same.runs: T or F whether speed and theta shoulf be the same on different trials

# N.cells <- 200; L  <- 200; L.total <- L; Ntrial <- 2; dt <- 0.001; N.x.basis <- 40; N.phi.basis <- 4; rate.ave <- 9; rate.bg = 0.1; seed <- 1; vx <- NULL; phi <- NULL; graphics <- T; adapt <- T; sample.template=F; randfields=F; randrates=F; fields.middle=F; var.trials=1/100000; same.runs=T; mu.v <- 6; i.clust <- 51:150

	set.seed(seed)

	sd.v <- 0.1 * var.trials # standard deviation of OU-speed (default: 0.1)
	sd.phi <- 2.5 * var.trials # standard deviation of theta-speed (default: 2.5)

	## generate a trajectory in 1D
	if (is.null(vx)){
		vx <- list(x=list(seq(0, by=0.6, length=435) %% 200), v=list(rep(600, 435)))
		for (i in 2:Ntrial) {
			vx$x[[i]] <- vx$x[[1]]
			vx$v[[i]] <- vx$v[[1]]
		}
	}

	tStartEnd <- matrix(NA, Ntrial, 2)
	times <- list()
	t0 <- dt
	for (i.trial in 1:Ntrial) {
		ttimes <- seq(t0, by=dt, length=length(vx$x[[i.trial]]))
		times[[i.trial]] <- ttimes
		tStartEnd[i.trial,] <- range(ttimes)
		t0 <- ceiling(max(ttimes + 1))
	}
	
	## generate ripple oscillation
	if (is.null(phi)){
		phi <- list()
		for (i.trial in 1:Ntrial){
			Tmax.i <- length(vx$x[[i.trial]]) * dt
			if (same.runs) sseed <- seed + 1 else sseed <- seed + 1 + i.trial
			pp <- sim.run.circ(Tmax = Tmax.i, xmax=2*pi, dt=dt, mu.v=150*2*pi, sd.v=sd.phi*2*pi, tau.v=0.5, seed=sseed, v.max=200*2*pi, v.min=5*2*pi, init.v=NULL)$x
			phi[[i.trial]] <- cbind(theta=NA*pp, amplitude=NA*pp, phase=pp)
		}
	}
	
	##################################################
	## 2. generate basis functions
	N.basis <- N.x.basis + N.phi.basis
	
	## 1D basis functions:
	## position, theta
	xx <- seq(0,L.total)
	x.basis <- gen.Gauss.basis.mat(x=xx, sds=5, n=N.x.basis, circular=T) # was 10 before...
	# if (graphics) {par(mfcol=c(2,1));	matplot(t(x.basis), t="l", col=rainbow(N.x.basis, end=0.7), lty=1, axes=F, xlab="distance (cm)", ylab="activation", main="spatial basis"); axis(1); axis(2, las=2)}

	L.phi.basis <- 25
	pphi <- seq(0,2*pi, length=L.phi.basis)
	phi.basis <- gen.basis.circ.mat(x=pphi, n=N.phi.basis, F)
	if (graphics) {matplot(1:25 / 25*360, t(phi.basis), t="l", col=rainbow(12), lty=1, axes=F, xlab="ripple phase (deg)", ylab="activation", main="ripple basis"); axis(1); axis(2, las=2)}

	## generate basis functions for connectivity / adaptation
	N.c.basis <- 6
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }
	
	xphi.basis <- matrix(NA, N.x.basis + N.phi.basis, ncol(x.basis)*4)
	xphi.basis[1:40,] <- x.basis
	xphi.basis[41:44,1:201] <- phi.basis[,1]
	xphi.basis[41:44,1:201+201] <- phi.basis[,7]
	xphi.basis[41:44,1:201+201*2] <- phi.basis[,13]
	xphi.basis[41:44,1:201+201*3] <- phi.basis[,19]
	
	#############################
	## generate weights for the basis
	
	N.allcells <- N.cells + length(i.clust)
	
	
	set.seed(seed)
	ww <- matrix(NA, N.allcells, N.x.basis + N.phi.basis)

	## generate weigths for the location -  basis functions
	w0 <- rep(NA, N.allcells)
	# w.template <- c(0, 2, 2,1.5,1.8,2, 2.5,2.5, 0,0)
	w.template <- c(1, 1.2, 1.5, 1.7, 1.8, 1.8, 1.7, 1.5,1.2,1)/4
	wfg.template <- c(1.6, 1.6, 1.7, 1.7, 1.8, 1.8, 1.7, 1.7,1.6,1.6)/4
	nx <- length(w.template)
	if (fields.middle) max.shift <- N.x.basis - nx else max.shift <- N.x.basis
	if (randfields) i.shifts <- sample(seq(0, max.shift-1), N.cells, T) else i.shifts <- (1:N.cells) %% max.shift
	i.shifts <- sort(i.shifts)

	if (is.null(phi.template)) phi.template <- 	c(0,0.4,1.1,0.4)

	for (i.cell in 1:N.cells){
		i.shift <- i.shifts[i.cell]

		if (i.cell %in% i.clust){ # if the cell is in the clustered group, we decrease its bg firing to almost zero, and only the place field firing remains
			rate.factor <- 0.38429
			w.mat <- rep(0, N.x.basis+nx)
			w.mat[(1:nx)+i.shift] <- wfg.template *4.7 - 1
			w.mat[1:nx] <- w.mat[1:nx] + w.mat[(1:nx)+N.x.basis]
			wLoc <- c(w.mat[1:N.x.basis], phi.template) - 1/2
			ratemap <- t(wLoc %*% xphi.basis)
			# points(exp(ratemap), col=3)

			
		} else { # bg + place field firing
			rate.factor <- 1
			w.mat <- rep(0, N.x.basis+nx)
			w.mat[(1:nx)+i.shift] <- w.template
			w.mat[1:nx] <- w.mat[1:nx] + w.mat[(1:nx)+N.x.basis]
			wLoc <- c(w.mat[1:N.x.basis], phi.template)
			ratemap <- t(wLoc %*% xphi.basis)
			# plot(exp(ratemap), ylim=c(0, 10))
		}

		if (randrates) r0 <- rgamma(1, 6*rate.ave, 6) * rate.factor else r0 <- rate.ave * rate.factor
		mean.rate <- mean(exp(ratemap))
		w0[i.cell] <- log(r0 / mean.rate)
		ww[i.cell,] <- wLoc
	}
	
	if (!is.null(i.clust)){ # background firing - no place field
		for (i.cell in (1+N.cells):N.allcells){
			rate.factor <- 0.61571
			w.mat <- rep(0, N.x.basis+nx)
			wLoc <- c(w.mat[1:N.x.basis], phi.template) - 0.05
			ratemap <- t(wLoc %*% xphi.basis)
			if (randrates) r0 <- rgamma(1, 6*rate.ave, 6) * rate.factor else r0 <- rate.ave * rate.factor
			mean.rate <- mean(exp(ratemap))
			w0[i.cell] <- log(r0 / mean.rate)
			ww[i.cell,] <- wLoc
		}
	}	
	# points(exp(ratemap.bg), col=2)
	# lines(exp(ratemap.fg) + exp(ratemap.bg), col=5, lwd=2)
	
	cat("w generated \n")
	# plot(exp(ratemap[1:200]))
	
	## basis activations and firing rate in the function of time
	raster <- list()
	spiketimes <- list()
	sp.cells <- rep(0, N.allcells)
	tt.cells <- 0
	if (graphics) rate.t <- rep(0, N.allcells)
	
	for (i.trial in 1:Ntrial){
		
		cat('starting trial',  i.trial, ' \n')
		cat('phi, ')
		phi.trial <- phi[[i.trial]][,3]
		cat('TT, ')
		TT <- length(vx$x[[i.trial]])
		cat('x.t with TT:', TT, '\n')
		x.t <- matrix(NA, N.basis, TT)

		for (t in 1:TT) {
			i.x <- which.min((xx - vx$x[[i.trial]][t])^2)
			i.phi <- which.min((pphi - phi.trial[t])^2)
			x.t[,t] <- c(x.basis[,i.x],phi.basis[,i.phi])
		}
		if (graphics) rate.t <- cbind(rate.t, exp(ww %*% x.t) * exp(w0))
		cat('basis functions calculated \n')
		
		# no connections
		N.wc <- 1
		w.c <- matrix(0, N.wc, N.c.basis)
		cons <- matrix(NA, 1, 2)
		
		# adaptation
		if (adapt){
			if (N.c.basis != 6) stop("w.self should be provided")
			w.burst <- c(-4, 1, 2, 0, -1, 0)
			w.adapt <- c(-4, 0, 0, -1/8, -1/4, -1/4)
			w.regular <- c(-4, -2, -1, -1, -1, 0)		
			w.refr <- c(-4, -1, 0, 0, 0, 0)
			w.self <- w.refr
			w.self <- matrix(rep(w.self, N.allcells), N.allcells, N.c.basis, byrow=T)	
		} else {
			w.self <- matrix(0, N.allcells, N.c.basis)
		}
		
		cat('starting the simulations for trial',  i.trial, ' \n')
		set.seed(seed * 1000 + i.trial)
		resp <- sim.glm(x.t, ww, cons, c.basis, w.c, w.self, w0, dt, rate.bg=rate.bg)
		cat("collecting the spikes \n")
		spt <- Matrix(0, N.allcells, TT)
		for (i in 1:nrow(resp$sp)) spt[resp$sp[i,1], resp$sp[i,2]] <- spt[resp$sp[i,1], resp$sp[i,2]] + 1
		cat('spike raster prepared... \n')
		raster[[i.trial]] <- spt
		spiketimes[[i.trial]] <- resp$sp[,c(2,1)]
		spiketimes[[i.trial]][,1] <- spiketimes[[i.trial]][,1] * dt + tStartEnd[i.trial,1]
		sp.cells <- sp.cells + rowSums(spt)
		tt.cells <- tt.cells + ncol(spt) * dt
		cat("\n trial", i.trial, " finished. \n")
	}
	
	print("simulations finished")
	if (graphics) rate.t <- rate.t[,-1]

	acells <- matrix(1, N.allcells, 4, dimnames=list(NULL, c("tetrode", "cell", "number of spikes", "firing rate")))
	acells[,2] <- 1:N.allcells
	acells[,3] <- sp.cells
	acells[,4] <- sp.cells / tt.cells
	
	
	syn_run <- list(runtype=1, t=tStartEnd, times=times, dists=vx$x, speed.1D=vx$v, raster=raster, acells=acells, lfp=list(phi), areas="CA1", spt=spiketimes)

	############################

	cat("graphics:", graphics, "\n")
	
	if (graphics){
		x <- unlist(vx$x)
		rasts <- raster[[1]]; for (i in 2:Ntrial) rasts <- cbind(rasts, raster[[i]])
		theta <- phi[[1]][,3]; for (i in 2:Ntrial) theta <- c(theta, phi[[i]][,3])
		sp <- raster2sp(rasts)
		Tmax <- length(x) * dt
		
		layout(matrix(c(1,4,2,3), 2), widths=c(1,4))
		par(mar=c(4,4,1,1))
		
		counts.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$counts * dt
		mids.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.allcells, 20)
		for (i in 1:N.allcells) rates.measured[i,] <- hist(x[rasts[i,] > 0], br=seq(0, max(x), length=21), plot=F)$counts / counts.x 
		matplot(mids.x, t(rates.measured), t="l", lty=2, col=rainbow(N.allcells, end=0.7), lwd=1, xlab="distance (cm)", ylab="rate (Hz)", axes=F); axis(1); axis(2, las=2)
		# matplot(xx[-1], t(ratemap), t="l", col=rainbow(N.allcells, end=0.7), lty=1, add=T)
		lines(mids.x, colMeans(rates.measured), t="l", lwd=2)
		abline(h=rate.ave, col=2)
		
		par(mar=c(1,3,1,4))
		kk <- nrow(spiketimes[[1]])
		plot(sp[1:kk,2], sp[1:kk,1], xlim=c(0, 400), pch="|", axes=F, xlab="", ylab=""); axis(2, las=2); mtext("cells", 2, 2, cex=0.85)
		lines(x*N.allcells/L.total, t="l", col=2); mtext("position (cm)", 4, col=2, line=2); axis(4, c(0,50, 100)*N.allcells/L.total, c(0,50,100), col=2)

		par(mar=c(4,3,1,4))
		t <- seq(dt, Tmax, by=dt)
		matplot(t[1:401], t(rate.t[, 1:401] + rate.bg), t="l", lty=1, col=rainbow(N.allcells, end=0.7), lwd=1, axes=F, xlab="time (s)", ylab=""); axis(2, las=2); axis(1); mtext("rate (Hz)", 2, 2, cex=0.85)
		
		counts.phi <- hist(theta, br=seq(0,2*pi, length=11), plot=F)$counts * dt
		mids.phi <- hist(theta, br=seq(0,2*pi, length=11), plot=F)$mids
		rates.measured <- matrix(NA, N.allcells, 10)
		for (i in 1:N.allcells) rates.measured[i,] <- hist(theta[rasts[i,] > 0], br=seq(0,2*pi, length=11), plot=F)$counts / counts.phi 
		par(mar=c(4,4,1,1))
		matplot(mids.phi, t(rates.measured), t="l", lty=2, col=rainbow(N.allcells, end=0.7), lwd=1, axes=F, xlab="ripple phase", ylab="rate (Hz)"); axis(1); axis(2, las=2)
		lines(mids.phi, colMeans(rates.measured), t="l", lwd=2)
			
	}
	syn_run

}


## square pulse stimulus + adaptation  - just to test adaptation kernels
gen.adaptcells <- function(Tmax =4, dt=0.001, graphics=T, N.cells=10, N.c.basis=6, rate.ave=10, w.self=NULL, rseed=22){
	# Tmax <- 2; dt <- 0.001; graphics <- T; rate.ave <- 10; N.cells <- 10; w.self <- NULL
	L <- Tmax / dt#; L1 <- round(0.15*L); L2 <- round(0.35*L); L3 <- round(0.65*L); L4 <- round(0.85*L)
	## stimulus
	x.t <- matrix(0, 1, L)
	# x.t[1, L1:L2] <- 1
	# x.t[1, L3:L4] <- 1
	w.x <- matrix(2, N.cells, 1)

	ratemap <- t(w.x %*% x.t)
	mean.rate <- mean(exp(ratemap))
	w0 <- log(rate.ave / mean.rate) 
	
	## generate basis functions for connectivity
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17
	}
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	if (graphics) {
		matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	}
	# no connections
	cons <- matrix(NA, 1, 2)
	N.wc <- 1
	w.c <- matrix(0, N.wc, N.c.basis)

	if (is.null(w.self)){
		if (N.c.basis != 6) stop("w.self should be provided")
		w.burst <- c(-4, 1, 2, 0, -1, 0)
		w.adapt <- c(-4, 0, 0, -1/8, -1/4, -1/4)
		w.regular <- c(-4, -2, -1, -1, -1, 0)
		w.refr <- c(-4, -1, 0, 0, 0, 0)
		w.self <- w.refr
		# w.self <- w.burst
	}
	
	w.self <- matrix(rep(w.self, N.cells), N.cells, N.c.basis, byrow=T)	

	if (graphics) {
		plot((w.self %*% c.basis)[1,], t="l", axes=F, ylab="gain", xlab="timestep")
		abline(h=0); axis(1); axis(2, las=2)
	}
	set.seed(rseed)
	resp <- sim.glm(x.t, w.x, cons, c.basis, w.c, w.self, w0, dt, rates=T)

	vm <- log(resp$r)
	for (ii in 1:nrow(resp$sp)) vm[resp$sp[ii,1], resp$sp[ii,2]] <- 10
	if (graphics) matplot(t(vm), t="l", lty=1)

	spiketimes <- resp$sp[,c(2,1)]
	spiketimes[,1] <- spiketimes[,1] * dt

	xs <- list(sp=spiketimes, x=x.t, r=resp$r, LL=resp$LL)	
	xs
}


##############################################################################
## location dependent cells with cell-to-cell connectivity, adaptation - no phase precession
gen.corrcells <- function(N.cells=10, L =100, Ntrial=10, dt=0.01, N.x.basis=11, rate.ave=2, seed=23, vx=NULL, phi=NULL, graphics=T, adapt=F, WW=NULL){
# N.cells <- 100; L <- 100; Ntrial <- 10; dt <- 0.001; seed = 52; N.x.basis <- 11; rate.ave <- 2; adapt <- F; WW <- WW; graphics <- T; vx <- NULL; phi <- NULL
	set.seed(seed)
	if (is.null(WW)) stop("connection matrix must be specified!") else {
		if (ncol(WW) != N.cells) stop("number of cells N.cells is not equal to the columns in the WW!")
		adapt <- T
	}
	
	## generate a trajectory in 1D
	if (is.null(vx)){
		vx <- sim.trial(N = Ntrial, xmax=L/100, dt =dt, seed=seed, sd.v = 0.1, mu.v=0.2, tau.v=1, v.min=0.1, v.max=0.5)
	}
	
	# generate time stamps for the trials
	tStartEnd <- matrix(NA, Ntrial, 2)
	times <- list()
	t0 <- dt
	for (i.trial in 1:Ntrial) {
		ttimes <- seq(t0, by=dt, length=length(vx$x[[i.trial]]))
		times[[i.trial]] <- ttimes
		tStartEnd[i.trial,] <- range(ttimes)
		t0 <- ceiling(max(ttimes + 1))
	}
	
	## generate theta oscillation
	if (is.null(phi)){
		phi <- list()
		for (i.trial in 1:Ntrial){
			Tmax.i <- length(vx$x[[i.trial]]) * dt
			pp <- sim.run.circ(Tmax = Tmax.i, xmax=2*pi, dt=dt, mu.v=8*2*pi, sd.v=2.5*2*pi, tau.v=0.5, seed=seed+1, v.max=12*2*pi, v.min=5*2*pi, init.v=NULL)$x
			phi[[i.trial]] <- cbind(theta=NA*pp, amplitude=NA*pp, phase=pp)
		}
	}
	
	## generate basis functions for place cells
	xx <- seq(0,L)
	x.basis <- gen.Gauss.basis.mat(x=xx, n=N.x.basis, sds=10)
	# if (graphics) matplot(t(x.basis), t="l", col=rainbow(12), lty=1, axes=F, xlab="distance (cm)", ylab="activation", main="spatial basis"); axis(1); axis(2, las=2)
	
	## generate basis functions for connectivity
	if (dt == 0.01) {
		min.phi <- 3
		max.phi <- 9
	}
	if (dt == 0.001) {
		min.phi <- 4
		max.phi <- 17.5
	}
	N.c.basis <- 6
	c.basis <- gen.cos.basis.mat(t=seq(1,0.25/dt), phis=seq(min.phi, max.phi, len=N.c.basis))
	# if (graphics) {
		# matplot(t(c.basis), t="l", col=rainbow(8), lty=1, axes=F, xlab="time (ms)", ylab="activation", main="coupling basis"); axis(1); axis(2, las=2)
	# }
		
	## generate random weigths for the basis functions
	w.xx <- matrix(rexp(N.cells*N.x.basis, 1), N.cells, N.x.basis)
	w.x <- sigm(w.xx, a=5, sl=5, th=2.3) # 90% around 0
	for (i in 1:N.cells) {
		if (max(w.x[i,]) < 2) {
			i.max <- which.max(w.x[i,])
			w.x[i, i.max] <- runif(1, 2,5)
		}
	}
	# hist(w.x)
	
	## weight of the constant basis to scale the firing rate 
	ratemap <- t(w.x %*% x.basis)
	i.w.x <- sort(apply(ratemap, 2, which.max), ind=T)$ix
	w.x <- w.x[i.w.x,]
	ratemap <- t(w.x %*% x.basis)
	
	mean.rate <- colMeans(exp(ratemap))
	w0 <- log(rate.ave / mean.rate) 

	ratemap <- t(exp(ratemap)) * exp(w0)
	# matplot(t(ratemap), t="l", lty=1, col=rainbow(120), lwd=2, axes=F, xlab="distance (cm)", ylab="rate (Hz)"); axis(1); axis(2, las=2)
	
	## basis activations and firing rate in the function of time
	raster <- list()	
	spiketimes <- list()
	sp.cells <- rep(0, N.cells)
	tt.cells <- 0
	if (graphics) rate.t <- rep(0, N.cells)
	for (i.trial in 1:Ntrial){
		TT <- length(vx$x[[i.trial]])
		x.t <- matrix(NA, N.x.basis, TT)
		for (t in 1:TT) {
			i.x <- which.min((xx - vx$x[[i.trial]][t])^2)
			x.t[,t] <- x.basis[,i.x]
		}
		if (graphics) 	rate.t <- cbind(rate.t, exp(w.x %*% x.t) * exp(w0))
		
		# connections
		N.wc <- sum(WW!=0)
		w.c <- matrix(0, N.wc, N.c.basis)
		w.exc <- c(0.2, 1, 0.2, -1/2, 0, 0)
		w.inh <- c(-2, -2, 0, 0, 0, 0)
		# plot(t(c.basis) %*% w.exc, t="l", lty=1); abline(v=c(0, 20, 40, 60), h=0)	
		# plot(t(c.basis) %*% w.inh, t="l", lty=1); abline(v=c(0, 20, 40, 60), h=0)	
		# matplot(t(c.basis), t="l", lty=1, col=rainbow(8))
		
		# pdf(file="cons.pdf", 5, 3)
		# plot(t(c.basis) %*% w.exc, t="l", lty=1, axes=F, log="", lwd=2, xlab="time (ms)", ylab="gain", col=2, ylim=c(-3,1)); axis(1); axis(2)
		# lines(t(c.basis) %*% w.inh, t="l", lty=1, col=4, lwd=2)
		# lines(t(c.basis) %*% w.adapt, t="l", lty=1, col=3, lwd=2); abline(h=0)
		# dev.off()
		
		i.con <- 1
		connections <- matrix(NA, N.wc, 2)
		for (i.pre in 1:N.cells){
			for (i.post in 1:N.cells){
				if (WW[i.post, i.pre] !=0){
					connections[i.con,] <- c(i.post, i.pre)
					if (WW[i.post, i.pre] > 0) w.c[i.con,] <- w.exc else w.c[i.con,] <- w.inh
					i.con <- i.con + 1
				}
			}
		}
				
		# adaptation
		if (adapt){
			if (N.c.basis != 6) stop("w.self should be provided")
			w.burst <- c(-6, 1, 2, 0, -1, -2)
			w.adapt <- c(-6, 0, 0, 0, -1/5, -1/5)
			w.regular <- c(-6, -2, -1, -1, -1, 0)		
			w.self <- w.adapt
			# plot(t(c.basis) %*% w.self, t="l", lty=1); abline(v=c(0, 20, 40, 60), h=0)	
			w.self <- matrix(rep(w.self, N.cells), N.cells, N.c.basis, byrow=T)	
		} else {
			w.self <- matrix(0, N.cells, N.c.basis)
		}

		# if (noadapt){
			# w.self <- matrix(0, N.cells, N.c.basis)
			# w.c <- matrix(0, N.wc, N.c.basis)		
		# }
				
		resp <- sim.glm(x.t, w.x, connections, c.basis, w.c, w.self, w0, dt)
		# plot(resp$sp[,2], resp$sp[,1], xlim=c(0, 1000), ylim=c(0, 20) , pch="|")
		# plot(resp$sp[,2], resp$sp[,1] , pch="|")
		
		spt <- matrix(0, N.cells, TT)
		for (i in 1:nrow(resp$sp)) spt[resp$sp[i,1], resp$sp[i,2]] <- spt[resp$sp[i,1], resp$sp[i,2]] + 1
		raster[[i.trial]] <- spt
		spiketimes[[i.trial]] <- resp$sp[,c(2,1)]
		spiketimes[[i.trial]][,1] <- spiketimes[[i.trial]][,1] * dt + tStartEnd[i.trial,1]
		sp.cells <- sp.cells + rowSums(spt)
		tt.cells <- tt.cells + ncol(spt) * dt
	}
	if (graphics) rate.t <- rate.t[,-1]

	acells <- matrix(1, N.cells, 4, dimnames=list(NULL, c("tetrode", "cell", "number of spikes", "firing rate")))
	acells[,2] <- 1:N.cells
	acells[,3] <- sp.cells
	acells[,4] <- sp.cells / tt.cells
	
	
	syn_run <- list(runtype=1, t=tStartEnd, times=times, dists=vx$x, speed.1D=vx$v, raster=raster, acells=acells, lfp=list(phi), areas="CA1", spt=spiketimes)

	if (graphics){
		x <- unlist(vx$x)
		rasts <- raster[[1]]; for (i in 2:Ntrial) rasts <- cbind(rasts, raster[[i]])
		theta <- phi[[1]][,3]; for (i in 2:Ntrial) theta <- c(theta, phi[[i]][,3])
		sp <- raster2sp(rasts)
		Tmax <- length(x) * dt
		
		layout(matrix(c(1,4,2,3), 2), widths=c(1,4))
		par(mar=c(4,4,1,1))
		
		print("matplot")
		matplot(t(ratemap), t="l", lty=1, col=rainbow(12), lwd=2, axes=F, xlab="distance (cm)", ylab="rate (Hz)"); axis(1); axis(2, las=2)
		counts.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$counts * dt
		mids.x <- hist(x, br=seq(0, max(x), length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(x[rasts[i,] > 0], br=seq(0, max(x), length=21), plot=F)$counts / counts.x 
		matplot(mids.x, t(rates.measured), t="l", lty=2, col=rainbow(12), lwd=1, add=T)
				
		par(mar=c(1,3,1,4))
		plot(sp[,2], sp[,1], pch="|", axes=F, xlab="", ylab=""); axis(2, las=2); mtext("cells", 2, 2, cex=0.85)
		lines(x/10, t="l", col=2); mtext("position (cm)", 4, col=2, line=2); axis(4, c(0,50, 100)/12+1, c(0,50,100), col=2)
		par(mar=c(4,3,1,4))
		t <- seq(dt, Tmax, by=dt)
		matplot(t, t(rate.t), t="l", lty=1, col=rainbow(12), lwd=2, axes=F, xlab="time (s)", ylab=""); axis(2, las=2); axis(1); mtext("rate (Hz)", 2, 2, cex=0.85)
		
		
		counts.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$counts * dt
		mids.phi <- hist(theta, br=seq(0,2*pi, length=21), plot=F)$mids
		rates.measured <- matrix(NA, N.cells, 20)
		for (i in 1:N.cells) rates.measured[i,] <- hist(theta[rasts[i,] > 0], br=seq(0,2*pi, length=21), plot=F)$counts / counts.phi 
		par(mar=c(4,4,1,1))
		matplot(mids.phi, t(rates.measured), t="l", lty=1, col=rainbow(12), lwd=1, axes=F, xlab="theta phase", ylab="rate (Hz)"); axis(1); axis(2, las=2)

	}
	
	syn_run
}


###################################################################
## generates 1-dimensional idealised place cells 
gen.wLoc <- function(X.basis, L.track, wmax=c(4, 5), alpha=0.57, beta=1/0.14, rate.ave=1, active.cells = F, seed=317){
	# input:
	# N.basis: number of place cells
	# L.track: length of the track in cm
	# wmax: max of w weights for the basis functions
	# alpha, beta: parameters of gamma distribution of place field propensity from Rich et al, 2014
	# rate.ave: average firing rate - 1 Hz
	# active.cells: simulate only cells that are active on this environment (at least 1 place fied)
	#N.basis <- 41; L.track <- 400; wmax <- c(4, 5); alpha <- 0.57; beta <- 1/0.14; active.cells <- T; seed <- 317
	######################################
	# output:
	# w: vector of weights for the basis functions
	N.basis <- nrow(X.basis)
	wLoc <- rep(0, N.basis)
	s.basis <- (L.track/100) / (N.basis-1)  # now in m 
	set.seed(seed)
	
	pf.propensity <- 2 # must be smaller than 1 pfs / m - Rich et al, 2014
	new.cell <- T
	while (new.cell){
		pf.propensity <- rgamma(1, alpha, beta) # place fields / m
		pf.propensity <- pf.propensity * s.basis # place fields / basis
		n.fields <- rpois(1, N.basis * pf.propensity)
		if (active.cells) {
			if (n.fields > 0)  new.cell <- F
		} else {
			new.cell <- F
		}
	}
	
	r0 <- rgamma(1, 6*rate.ave, 6)
	r0 <- n.fields * r0  + 1/10 # average firing rate, Hz

	if (n.fields > 0){
		i.ws <- sample(seq(1,N.basis), n.fields)
		for (i.w in i.ws){
			wLoc[i.w] <- runif(1, wmax[1], wmax[2])
			if (i.w > 1) {if (runif(1)>1/2) wLoc[i.w-1] <- runif(1, wmax[1], wmax[2])}#wLoc[i.w] * runif(1, 0, 1)}
			if (i.w < N.basis) if (runif(1)>1/2) {wLoc[i.w+1] <- runif(1, wmax[1], wmax[2])}#wLoc[i.w] * runif(1, 0, 1)}
		}
	}
	
	ratemap <- t(wLoc %*% X.basis)
	mean.rate <- mean(exp(ratemap))
	w0 <- log(r0 / mean.rate) 
	ratemap <- exp(ratemap) * exp(w0)
	attributes(wLoc) <- list(w0=w0)
	wLoc	
}

