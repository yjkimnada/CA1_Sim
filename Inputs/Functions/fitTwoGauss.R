## we fit the orientation tuned responses with two Gaussians + a baseline
## these are used to derive parameters of the cell's responses 

dcGauss <- function(x, mu, sigma){
	# poor man's version of the circular normal - not normalised, but the peak is 1
	y <- dnorm(x, mu, sigma) + dnorm(x, mu + 2*pi, sigma)+ dnorm(x, mu - 2*pi, sigma)
	factor <- dnorm(mu, mu, sigma) + dnorm(mu, mu + 2*pi, sigma)+ dnorm(mu, mu - 2*pi, sigma)
	y <- y / factor
	y
}

# plot(lambda, dcGauss(lambda, pi, 0.36))

twoGauss <- function(pars, x){
	# baseline + two circular Gauss
	sigma = exp(pars[1])
	a1 = exp(pars[2]) # amplitude at mu
	a2 = exp(pars[3]) # amplitude at mu + pi
	b = exp(pars[4]) # baseline
	mu <- pars[5]
	if (mu > pi) mu2 <- mu - pi else mu2 <- mu + pi
	y <- b + a1 * dcGauss(x, mu, sigma) + a2 * dcGauss(x, mu2, sigma)
	y
}
# pars <- c(-1, 3, 2, 1, pi*3/2)
# plot(lambda, twoGauss(pars, lambda), t='o')

errTwoGauss <- function(pars, lambda, resp, alpha=0){	
	resp.est <- twoGauss(pars, lambda)
	err <- mean((resp-resp.est)^2) + alpha * pars[1]^2
	err	
}


#####################################

oneGauss <- function(pars, x){
	# baseline + two circular Gauss
	sigma = exp(pars[1])
	a1 = exp(pars[2]) # amplitude at mu
	b = exp(pars[3]) # baseline
	mu <- pars[4]
	if (mu < 0) mu <- mu %% (2*pi)
	if (mu > (2*pi)) mu <- mu %% (2*pi)
	y <- b + a1 * dcGauss(x, mu, sigma)
	y
}
# pars <- c(-1, 3, 2, 1, pi*3/2)
# plot(lambda, twoGauss(pars, lambda), t='o')

errOneGauss <- function(pars, lambda, resp, alpha=0){	
	resp.est <- oneGauss(pars, lambda)
	err <- mean((resp-resp.est)^2) + alpha * pars[1]^2
	err	
}
