### functions to generate synthetic population activity in the hippocampus
### firing of place cells while the rat explores a N-Dimensional environment
### 1. sim.run: motion of the animal in N-dimensions
### 2. check.boundary: check whether the boundary of the n-dimensional hyper-cube has been reached. If yes, the animal bounces back from the wall.


sim.run <- function(Tmax=10, xmax=1, dt=0.01, ndim=2, mu.v=0.4, sd.v=0.1, tau.v=0.1, m.dv=0.02, seed=123, v.max=Inf, v.min=-Inf, init.v=NULL){

# # Simulating the movement of the rat in an n-dimensional box
# # the speed is an OU process, smoothly changing with a mean of mu.v
# # the direction is changed in each timestep by adding a random vector to the current velocity 
# # and then renormalising its length
#######################################
# # input:
# # Tmax: duration of the simulations (s)
# # xmax: the size of the n-dimensional square box (m)
# # dt: time resolution (s)
# # ndim: the dimensionality of the random walk
# # mu.v: mean speed (m/s)
# # sd.v: sd of the speed 
# # tau.v: smoothness of the speed profile (s)
# # m.dv: the mean change of velocity
######################################
# # output: a list with elements: 
# # x: position (m)
# # v: velocity (m/s)

  	set.seed(seed)
	var <- sd.v^2
	
	if (mu.v + sd.v > v.max)  stop("mean + sd > v.max")
	if (mu.v - sd.v < v.min)  stop("mean - sd < v.min")
	
	# init the speed
	if (is.null(init.v)) speed.last <- rnorm(1, mu.v, sd.v) else speed.last <- init.v
	if (speed.last > v.max)  speed.last <- mu.v
	if (speed.last < v.min)  speed.last <- mu.v
	  
  	#### The hidden state variables
	if (length(xmax) > 1) x.last <- xmax/2 else x.last <- rep(xmax/2, ndim)# we start from the middle
	x.mat <- matrix(NA, ndim, Tmax/dt)
	x.mat[,1] <- x.last
	velocity.last <- rep(sqrt(speed.last^2/ndim), ndim) # velocity of leaving

	velocity.mat <- matrix(NA, ndim, Tmax/dt)
	velocity.mat[,1] <- velocity.last
		  	
	Q <- sqrt(var * 2 / tau.v)
 	for (i in 1:(Tmax/dt-1)){
		# 1. arriving to the new place
		# print(vv)
	    x.new <- x.last + velocity.last * dt
		vx <- check.boundary(velocity.last, x.new, xmax)
		velocity.last <- vx$v # this does not change the speed only the direction!
		x.new <- vx$x
		x.mat[,i+1] <- x.new

		# change the speed
    		speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v + Q * rnorm(1,0,1) * sqrt(dt)
    		if (speed.new < v.min)  speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v
	    if (speed.new > v.max)  speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v

		# change the velocity
		vv.add <- runif(ndim, -1,1)
		vv.add <- m.dv * vv.add / sqrt(sum(vv.add^2))

		velocity.new <- velocity.last + vv.add
		velocity.new <- velocity.new * speed.new / sqrt(sum(velocity.new^2)) # normalise to speed.new

		velocity.last <- velocity.new
		speed.last <- speed.new
		velocity.mat[,i+1] <- velocity.last
		x.last <- x.new
  	}
  	
  	vx <- list(x=x.mat, v=velocity.mat)
  	vx
}

########################################################################################

sim.trial <- function(N=10, xmax=1, dt=0.01, mu.v=0.4, sd.v=0.1, tau.v=0.1, m.dv=0.02, seed=123, v.max=1, v.min=0.05, init.v=NULL, out.cm=F){

# # Simulating the runs of the rat on a 1-dimensional box
# # the speed is an OU process, smoothly changing with a mean of mu.v
#######################################
# # input:
# # N: number of trials
# # xmax: the size of the n-dimensional square box (m)
# # dt: time resolution (s)
# # mu.v: mean speed (m/s)
# # sd.v: sd of the speed 
# # tau.v: smoothness of the speed profile (s)
# # m.dv: the mean change of velocity
# out.cm: output in centimeters
######################################
# # output: a list with elements: 
# # x: list of positions at each trial (m)
# # v: list of velocities at each trial (m/s)
# # N <- 10; xmax <- 1; dt <- 0.01; mu.v=0.2; sd.v=0.1; tau.v=1; m.dv=0.02; seed=123; v.max=0.6; v.min=0.05; init.v=NULL
  	set.seed(seed)
	var <- sd.v^2
	
	if (mu.v + sd.v > v.max)  stop("mean + sd > v.max")
	if (mu.v - sd.v < v.min)  stop("mean - sd < v.min")
	
	# init the speed
	xx <- list()
	vv <- list()
	
	for (i.trial in 1:N){
		if (is.null(init.v)) speed <- rnorm(1, mu.v, sd.v) else speed <- init.v
		if (speed > v.max)  speed <- mu.v
		if (speed < v.min)  speed <- mu.v
		  
	  	#### The hidden state variables
		x.last <- 0
		x.mat <- rep(NA, xmax / v.min / dt)
		x.mat[1] <- x.last
		speed.mat <- rep(NA, xmax / v.min / dt)
		speed.mat[1] <- speed
			  	
		Q <- sqrt(var * 2 / tau.v)
		
		i <- 1
		while (x.last < xmax){
			# 1. arriving to the new place
			# print(vv)
		    x.new <- x.last + speed * dt
			x.mat[i+1] <- x.new
	
			# change the speed
	    	speed <- speed + (mu.v - speed) * dt / tau.v + Q * rnorm(1,0,1) * sqrt(dt)
	    	while (speed < v.min)  speed <- speed + (mu.v - speed) * dt / tau.v
		    while (speed > v.max)  speed <- speed + (mu.v - speed) * dt / tau.v
	
			speed.mat[i+1] <- speed
			x.last <- x.new
			i <- i + 1
	  	}
	  	
	  	x.mat <- x.mat[1:i]
	  	speed.mat <- speed.mat[1:i]
	  	if (out.cm) scale.cm <- 100 else scale.cm <- 1
	  	xx[[i.trial]] <- x.mat  * scale.cm # change to cm from m
	  	vv[[i.trial]] <- speed.mat  * scale.cm # change to cm/s from m/s
	}  	
  	vx <- list(x=xx, v=vv)
  	vx
}

### check.boundary: check whether the boundary of the n-dimensional hyper-cube has been reached. If yes, the animal bounces back from the wall.

check.boundary <- function(v, x, xmax){
# # input: speed, position, boundaries and the time of the last turn
# # output: new speed vector
	i.turn <- which(((x>xmax) | (x < 0)))
	if (length(i.turn > 0)){
		for (i in i.turn){
			v[i] <- (-1) * v[i]
			if (x[i]>xmax[i]) x[i] <- 2 * xmax[i] - x[i] else x[i] <- (-1) * x[i]
		}
	}
	vx <- list(v=v, x=x)
	vx
}


# vv <- get.vv(10, c(pi/3, pi))
# sqrt(sum(vv^2)) - should be 10

# vx <- sim.run(Tmax = 20, xmax=c(5,1,1), ndim=3, seed=20, m.dv=0.05)
# plot(vx$x[1,], vx$x[2,], t="l")
# plot(vx$x[2,], vx$x[3,], t="l")

# matplot(t(vx$x), t="l", lty=1, col=rainbow(7))

# matplot(t(abs(vx$v)), t="l", lty=1)

# library(scatterplot3d)
# scatterplot3d(vx$x[1,],vx$x[2,],vx$x[3,], pch=".", color=rainbow(2000))

# L.track <- 12 # m
# dt <- 0.02 # 20 ms
# t <- seq(0, Tmax, by=dt)

# mu.v <- 0.40 # running speed - 40 cm/s
# sd.v <- 0.10 # m/s
# v.min <- 0.10 # m/s
# v.max <- 1 # m/s

# vx <- sim.run(L.track, dt, mu.v, sd.v, 0.5, seed=1, v.max, v.min)

# plot(vx$v, t="l")
# plot(vx$x, t="l")

# c(mean(vx$v), sd(vx$v))
# acf(vx$v, lag.max=200)

####################################################################
### sim.run.circ: simulates motion in a 1D circular track
### motion is smooth, as it is defined as an OU process in the speed

sim.run.circ <- function(Tmax=10, xmax=1, dt=0.01, mu.v=0.4, sd.v=0.1, tau.v=0.1, seed=123, v.max=Inf, v.min=-Inf, init.v=NULL){

# # Simulating the movement of the rat in a 1-dimensional bircle
# # the speed is an OU process, smoothly changing with a mean of mu.v
# # the direction is changed in each timestep by adding a random vector to the current velocity 
# # and then renormalising its length
#######################################
# # input:
# # Tmax: duration of the simulations (s)
# # xmax: the length of the 1-dimensional square circle (m)
# # dt: time resolution (s)
# # mu.v: mean speed (m/s)
# # sd.v: sd of the speed 
# # tau.v: smoothness of the speed profile (s)
######################################
# # output: a list with elements: 
# # x: position (m)
# # v: velocity (m/s)

  	set.seed(seed)
	var <- sd.v^2
	
	if (mu.v + sd.v > v.max)  stop("mean + sd > v.max")
	if (mu.v - sd.v < v.min)  stop("mean - sd < v.min")
	
	# init the speed
	if (is.null(init.v)) speed.last <- rnorm(1, mu.v, sd.v) else speed.last <- init.v
	if (speed.last > v.max)  speed.last <- mu.v
	if (speed.last < v.min)  speed.last <- mu.v
	  
  	#### The hidden state variables
	x.last <- xmax / 2 # we start from the middle
  	x.vec <- rep(NA, Tmax/dt)
  	x.vec[1] <- x.last
  	speed.vec <- rep(NA, Tmax/dt)
  	speed.vec[1] <- speed.last
		  	
	Q <- sqrt(var * 2 / tau.v)
 	for (i in 1:(Tmax/dt-1)){
		# 1. arriving to the new place
		# print(vv)
	    x.new <- x.last + speed.last * dt
		if (x.new > xmax) x.new <- x.new - xmax
	  	x.vec[i+1] <- x.new

		# change the speed
	    	speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v + Q * rnorm(1,0,1) * sqrt(dt)
   	 	if (speed.new < v.min)  speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v
	    if (speed.new > v.max)  speed.new <- speed.last + (mu.v - speed.last) * dt / tau.v

		# change the velocity
		speed.last <- speed.new
	  	speed.vec[i+1] <- speed.last
		x.last <- x.new
  	}
  	
  	vx <- list(x=x.vec, v=speed.vec)
  	vx
}


# vx <- sim.run.circ(Tmax = 20, xmax=2, seed=20)
# plot(vx$x, vx$v, t="l")
# plot(vx$x, t="l")

# vx <- sim.run.circ(Tmax = 20, xmax=2, seed=20, sd.v=0.25, v.min=0.1)
# plot(vx$x, vx$v, t="l")
# plot(vx$x, t="l")
