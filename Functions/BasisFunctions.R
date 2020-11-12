### defining basis functions
### basis functions are defined in two different ways: 
# # 1. by their parameters - separate function is needed to evaluate them 
# # 2. by a matrix providing their activity as the function of the variable (.mat functions)

#######################################################################
# generating a set of Gaussian bump basis functions in n-dimension
gen.Gauss.basis <- function(ranges=list(c(0, 200)), ns=c(11), force.sd=NULL, circular=NULL){
	# ranges: a list describing the ranges to be covered with the basis functions in each dimensions
	# ns: the number of basis function to be created in each dimension
	# the total number of basis functions will be prod(ns)

	ndim <- length(ns)
	if (is.null(circular)) circular <- rep(F, ndim)
	n.basis <- prod(ns)
	basis <- list(mu=list(), sd=list())
	sd <- unlist(lapply(ranges, max)) / ns # same sd for all basis

	if (!is.null(force.sd)){
		i.sd <- which(!is.na(force.sd))
		sd[i.sd] <- force.sd[i.sd]
	}

	mus <- list()
	for (idim in 1:ndim){
		if (circular[idim]) {
			dmus <- (ranges[[idim]][2] - ranges[[idim]][1]) / ns[idim]
		} else {
			dmus <- (ranges[[idim]][2] - ranges[[idim]][1]) / (ns[idim] - 1)
		}
		mus[[idim]] <- seq(ranges[[idim]][1], by=dmus, length=ns[idim])
	}
	mu <- expand.grid(mus)
	basis <- list(mu=mu, sd=sd)
	basis
}


eval.Gauss.basis <- function(x, basis, circular=F, period=NA){
	# x: the input at a given time - can be a vector
	# basis: list of the basis function parameters, 
	# basis$mu: k-dimensional mean
	# basis$sd: k-dimensional sd
	# circular: any of the dimensions should be treated as circular?
	# period: vector: NA for non-circular, otherwise the period of the circular dimension
	# output: the activation of the n basis at location x
	
	n <- nrow(basis$mu)
	f.x <- rep(NA, n)
	if (circular){
			j.circular <- which(!is.na(period))
			j.not.circular <- which(is.na(period))
		for (i in 1:n){
			f.xi.circular <- rep(NA, length(j.circular))
			k <- 1
			for (j in j.circular){
				f.xi.circular[k] <- exp(-1/2*(x[j]-basis$mu[i,j])^2/(basis$sd[j]^2)) + exp(-1/2*(x[j]-basis$mu[i,j]-period[j])^2/(basis$sd[j]^2)) + exp(-1/2*(x[j]-basis$mu[i,j]+period[j])^2/(basis$sd[j]^2))
				k <- k + 1
			}
			if (length(j.not.circular) > 0) f.x[i] <- prod(f.xi.circular) * exp(-1/2*sum( (x[j.not.circular]-basis$mu[i,j.not.circular])^2/(basis$sd[j.not.circular]^2))) else f.x[i] <- prod(f.xi.circular)
		}
	} else {
		for (i in 1:n){
			f.x[i] <- exp(-1/2*sum( (x-basis$mu[i,])^2/(basis$sd^2)))		
		}
	}
	f.x	
}

# gen.Gauss.basis(circular=T)

# ranges <- list(c(0, 100), c(0, 2*pi)); ns <- c(10, 2)
# basis <- gen.Gauss.basis(ranges, ns, circular=c(F,T), force.sd=c(NA, pi/8))

# xx <- seq(0,100)
# pphi <- seq(0, 2*pi, length=25)
# rates <- array(NA, dim=c(20, 25, 101))
# for (x in 1:101){
	# for (phi in 1:25){
		# rates[,phi, x] <- eval.Gauss.basis(c(xx[x], pphi[phi]), basis, circular=T, period=c(NA, 2*pi))
	# }
# }
# for (i.basis in 1:20){
	# image(t(rates[i.basis,,]))
	# readline(i.basis)
# }


# ranges <- list(c(0, 100)); ns <- c(20)
# basis <- gen.Gauss.basis(ranges, ns, force.sd=5)
# xx <- seq(0,100)
# rates <- array(NA, dim=c(ns, 101))
# for (x in 1:101){
	# rates[, x] <- eval.Gauss.basis(xx[x], basis)
# }
# matplot(t(rates), t="l", col=rainbow(ns), lwd=2, lty=1)


# matplot(t(rates[,13,]), t="l")
# matplot(t(rates[,,1]), t="l")


# matplot(t(gen.Gauss.basis()), t="l", col=rainbow(50), lty=1)
# matplot(t(gen.Gauss.basis(n=13)), t="l", col=rainbow(16), lty=1)


#######################################################################
# generating a set of Gaussian bump basis functions in 1-dimension - matrix version
gen.Gauss.basis.mat <- function(x=seq(0, 200), sds=5, n=(max(x)/sds) + 1, circular=F){
	# t: space in cm or time in s
	# sds: the sd of the bases
	# n: number of the basis functions to be generated

	mat.f <- matrix(0, n, length(x))

	if (circular) {
		mus <- seq(min(x), max(x), length=n+1)
		period <- diff(range(mus))
		mus <- mus[1:n]
	} else mus <- seq(min(x), max(x), length=n)
	
	
	
	for (i in 1:n){
		mat.f[i,] <- exp(-(x-mus[i])^2/(2*sds^2))
	}

	if (circular) {
		for (i in 1:n){
			mat.f[i,] <- mat.f[i,] + exp(-(x-mus[i] - period)^2/(2*sds^2))
			mat.f[i,] <- mat.f[i,] + exp(-(x-mus[i] + period)^2/(2*sds^2))
		}
	}
	mat.f
}

# matplot(t(gen.Gauss.basis.mat()), t="l", col=rainbow(50), lty=1)
# matplot(t(gen.Gauss.basis.mat(n=13)), t="l", col=rainbow(16), lty=1)


#######################################################################
## generates circular basis functions - phase preference of neurons
gen.basis.circ.mat <- function(x, n, graphics = F){
	# input:
	# x: input phases
	# n: number of basis
	# graphics: plot the basis
	######################################
	# output:
	# bb - basis activation in the function of phase
	
	if (min(x) != 0) x <- x - min(x)
	if (max(x) != (2*pi)) x <- x * 2 * pi / max(x)
	
	bb <- matrix(0, n, length(x))
	mus <- seq(0, by=2*pi/n, length=n)
	sd <- 2 * pi / n / 2

	for (ii in 1:n){
		xx <- 1 * exp(-(x-mus[ii])^2/(2*sd^2))
		xx <- xx + 1 * exp(-(x-mus[ii] + 2*pi)^2/(2*sd^2))
		xx <- xx + 1 * exp(-(x-mus[ii] - 2*pi)^2/(2*sd^2))
		bb[ii,] <- xx
	}

	if (graphics) matplot(t(bb), t="l")
	bb
}

# generating a set of raised cosine bump basis functions - as in Pillow et al., 2008
# b_j(t) = (1/2) * cos( a * log[t + c] + phi_j) + (1/2)
# t <- seq(0, 100, len=1000)
# a <- 3.75
# c <- 0.01
# M <- 10
# phis <- seq(-2,14, len=M)

#######################################################################
## cosyne basis for connectivity
gen.cos.basis.mat <- function(t=seq(1, 1000, len=1000), a=3.75, c=0.01, phis=seq(3,22, len=10), monoton=F, first=T){
	# t: time in ms
	# a: frequency of the basis
	# c: offset of the log
	# phis: phase of the cosine
	# monoton: monoton decreasing basis functions
	# first: the first basis is monoton decreasing
	M <- length(phis)
	mat.f <- matrix(0, length(phis), length(t))
	for (i in 1:length(phis)){
		phi <- phis[i]
		x <- a * log(t + c)
		f <- 0.5 * cos(x-phi) + 0.5
		f[x>phi+pi] <- 0
		f[x<phi-pi] <- 0
		if ((first) & (i==1)) {
			i.max <- which.max(f)
			f[1:i.max] <- max(f)
		}
		if (monoton){
			if (i>1){
				i.max <- which.max(f)
				f[1:i.max] <- max(f)			
			}
		}
		mat.f[i,] <- f		
	}
	mat.f
}

# # 
# matplot(t(gen.cos.basis.mat(t=seq(1,100), phis=seq(3,14,len=10))), t="l", col=rainbow(10), lty=1)
# matplot(t(gen.cos.basis.mat(t=seq(1,100), phis=seq(3,14,len=5))), t="l", col=rainbow(5), lty=1)
# matplot(t(gen.cos.basis.mat(t=seq(1,100), phis=seq(3,14,len=20))), t="l", col=rainbow(20), lty=1)


# definition of the logistic sigmoid function - useful for lot's of stuff!
sigm <- function(x, a=1, sl=1, th=0, d=0){
	# a / (1 + exp(-sl(x-th)))
	#a - amplitude: limit in inf; sl - slope; th - threshold
	# d: the derivative of the sigmoid function at x. d=0: the function; d=1: first derivative; d=2: second derivative
	# derivative with respect to the parameters sl or th
	bc <- exp(-sl*(x-th))
	if (is.numeric(d)){		
		if (d>2) return ("d must be smaller than 2!")
		if (d==0) {
			sx <- a / (1+bc)
			} else {
				if (d==1) {
					sx <- (a*sl*bc)/(1+bc)^2
				} else {
					sx <- (a*sl^2*bc*(bc-1))/(1+bc)^3
				}
				sx[(1+bc == bc)] <- 0
			}
		} else {
			if (d == "a"){
				sx <- 1 / (1+bc)
			}
			if (d == "b"){
				sx <- a * (x-th) * bc / (1+bc)^2	
			}
			if (d == "c"){
				sx <- -a * sl * bc / (1+bc)^2	
			}
			if (!(d %in% c("a", "sl", "th"))) return("invalid d value")
		}
	sx
}

########################################################################################
gen.W.basis <- function(alldist, allarm, sds=10, graphics = T, dir=F){
	# basis functions for the W-maze
	# allarm: vector of arm index -1, 0 or 1
	# alldist: vector of the 1D distances
	# sds: the sd of the bases
	# dir: place fields can be directional - this is implemented by having directional basis functions
	# sds <- 10; graphics <- T

	# plot(allpos[,1], alldist, pch=".")
	
	## generate basis functions for the total distance
	Xmax <- ceiling(max(alldist))
	Xmin <- floor(min(alldist))
	x <- seq(Xmin,Xmax)
	mat1 <- gen.Gauss.basis.mat(x=x, sds)
	n1 <- nrow(mat1)
	# matplot(x, t(mat1), col=viridis(n1, option="D"), t="l", lty=1, lwd=2)

	## find the distance of the junction
	d.junct <- floor(min(alldist[allarm==-1]))
	# abline(v=d.junct)
	i.junct <- which.min((x - d.junct)^2)

	## basis functions with peak after the junction are duplicated
	i.junct.basis <- which.max(mat1[,i.junct])+1
	mat2 <- mat1[(i.junct.basis:n1),]
	n2 <- nrow(mat2)
	# matplot(x, t(mat1[1:(i.junct.basis-1), ]), col=viridis(n1, option="D")[1:i.junct.basis], t="l", lty=1, lwd=2, ylim=c(-1,2))
	# abline(v=d.junct)
	# matplot(x, t(mat2)+1, col=viridis(n1, option="D")[i.junct.basis:n1], t="l", lty=1, lwd=2, add=T)
	# matplot(x, t(mat1[i.junct.basis:n1, ])-1, col=viridis(n1, option="D")[i.junct.basis:n1], t="l", lty=1, lwd=2, add=T)
	
	mat <- rbind(mat1[i.junct.basis:n1, ], mat1[1:(i.junct.basis-1), ], mat2)
	basis.arm <- rep(0, n1+n2)
	basis.arm[(n2+1):n1] <- 0
	basis.arm[1:n2] <- 1
	basis.arm[(n1+1):(n1+n2)] <- (-1)

	if (dir){
		mat <- rbind(mat, mat)
		basis.arm <- c(basis.arm, basis.arm)
		dir.arm <- rep(-1, length(basis.arm))
		dir.arm[1:(length(basis.arm)/2)] <- 1
	}

	if (dir) {
		attributes(mat) <- list(dim=dim(mat), arm=basis.arm, dir=dir.arm, x=x)
	} else {
		attributes(mat) <- list(dim=dim(mat), arm=basis.arm, x=x)
	}


	if (graphics){
		# print(dim(mat))
		# print(length(x))
		if (dir) {
			matplot(x, t(mat + seq(1,2*(n1+n2))/2), t="l", col=rainbow(8)[basis.arm+2 + 2*(dir.arm+1)], lty=1) 
		} else {
			matplot(x, t(mat + seq(1,(n1+n2))/2), t="l", col=rainbow(4)[basis.arm+2], lty=1) 
		}
		abline(v=d.junct)
	}
	
	mat

}



eval.W.basis <- function(dist.t, arm.t, mat){
	# arm.t: vector of arm index -1, 0 or 1
	# dist.t: vector of the 1D distances
	# dist.t <- alldist; arm.t <- allarm; mat <- W.basis
	x <- attributes(W.basis)$x
	arms <- attributes(W.basis)$arm
	dir.bases <- F
	if ("dir" %in% names(attributes(mat))) {
		dir.bases <- T
		dir.mat <- attributes(W.basis)$dir # directionality of the basis; -1 or 1
		dir.run <- sign(diff(dist.t)) # direction of the run; -1 or 1
		dir.run <- c(dir.run[1], dir.run)
	}

	act.W <- matrix(NA, length(dist.t), nrow(W.basis))	
	for (t in 1:length(dist.t)){
		dd <- dist.t[t]
		aa <- arm.t[t] # the arm - only basis with arm=0 and the current arm will be kept active!
		i.dd <- which.min((dd-x)^2)
		act.t <- W.basis[,i.dd] # all basis activated
		i.off.aa <- (arms != 0) & (arms != aa) # central arm activation is never switched off
		act.t[i.off.aa] <- 0
		if (dir.bases){
			i.off.dd <- (dir.mat != dir.run[t])
			act.t[i.off.dd] <- 0
		}
		act.W[t,] <- act.t
	}
	act.W
}
