## A script to simulate response of a population of neurons to external input with coupling terms and adaptation

sim.glm <- function(x.t, w.x, cons, c.basis, w.c, w.self, w0, dt=0.001, graphics=F, rates=F, rate.bg=0){
	# x.t: external input (K x T)
	# w.x: the tuning of the cells to different inputs, N x K
	# cons: the connectivity between cells; N x 2, post, pre
	# c.basis: basis functions for coupling
	# w.c: coupling weights
	# w.self: self-coupling or adaptation weights
	# w0: scaling the firing rate
	# dt: time resolutions for spiking
	# cons <- connections; graphics <- T; rates=F
	# rate.bg: background firing - mostly should be 0! Resistent to adaptation and bursting.
	
	TT <- ncol(x.t)
	T0 <- ncol(c.basis)
	Tmax <- T0 + TT
	N.cells <- nrow(w.x)
	spikes <- matrix(0, N.cells, Tmax)
	m.inp <- matrix(0, N.cells, 3); colnames(m.inp) <- c("input", "self", "peers")
	v.inp <- matrix(0, N.cells, 3); colnames(v.inp) <- c("input", "self", "peers")
	if (rates) rr <- matrix(NA, N.cells, Tmax)
	LL <- rep(0, N.cells) # log likelihood for the cells
	
	for (t in (T0 +1):Tmax){
		# the effect of spatial input
		# t <- t + 1
		inp.x <- w.x %*% x.t[,t-T0]
		m.inp[,1] <- m.inp[,1] + inp.x
		v.inp[,1] <- v.inp[,1] + inp.x^2
		
		# the effect of spatial past spikes
		past.sp <- spikes[,(t-1):(t-T0)]
		act.bases <- past.sp %*% t(c.basis)
		inp.peers <- rep(0, N.cells)
		inp.self <- rep(0, N.cells)
		
		for (i.post in 1:N.cells){
			inp.self[i.post] <- act.bases[i.post,] %*% w.self[i.post,]

			i.cons <- which(cons[,1] == i.post) # connections where i is postsynaptic
			if (length(i.cons) > 0 ){
				for (i.con in i.cons){
					i.pre <- cons[i.con,2]
					inp.peers[i.post] <- inp.peers[i.post] + act.bases[i.pre,] %*% w.c[i.con,]
				}
			}
		}

		m.inp[,2] <- m.inp[,2] + inp.self
		v.inp[,2] <- v.inp[,2] + inp.self^2
		m.inp[,3] <- m.inp[,3] + inp.peers
		v.inp[,3] <- v.inp[,3] + inp.peers^2

		lambda <- dt * exp(inp.peers) * exp(inp.self) * exp(inp.x) * exp(w0) + rate.bg * dt
		if (rates) rr[,t] <- lambda / dt
		p.silence <- exp(-lambda) # from Poisson
		sp <- runif(N.cells)
		sp[sp<p.silence] <- 0
		sp[sp>0] <- 1 # P(s=1) = P(s>0)!
		spikes[,t] <- sp
		if ((t %% 1000) == 0) cat(t, " ")

		# LL <- LL + dpois(sp, lambda, log=T)	
		LL <- LL + dbinom(sp, 1, 1- exp(-lambda), log=T)
	}
	sp <- cbind((which(!!spikes) - 1) %% N.cells + 1, (which(!!spikes) -1) %/% N.cells + 1)
	sp[,2] <- sp[,2] - T0
	if (graphics) {plot(sp[,2], sp[,1], pch="|", axes=F, xlab="time (ms)", ylab="cells"); axis(1)}
	m.inp <- m.inp/TT
	v.inp <- v.inp/TT - m.inp^2
	
	resp <- list(sp=sp, m.inp=m.inp, v.inp=v.inp, LL=LL)
	if (rates) {
		rr <- rr[,-(1:T0)]
		resp$r <- rr
	}
	resp		
}

