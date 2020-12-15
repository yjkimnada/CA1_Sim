########################################################
comb.lin <- function(resp, delay=NULL, n=NULL, L=length(resp), dt=0.1, times=NULL, ampls=NULL){
## this function adds its resp arguments n times with delay (delay, n). It returns the first L point, starting at 0.
## alterntively you can specify the exact spacement of the functions to add up (times) - see second part

	resp <- resp - resp[1]
	if (is.null(times)){
		if (is.null(delay)) stop("delay or times must be specified")
		if (is.null(n)) stop("n or times must be specified")
		print(delay)
		delay <- delay/dt
		l <- length(resp)
		if (L > l) stop ("L is too large")
		if (delay * (n-1) >= l) stop ("delay is too large")
		if (n == 1) return(resp)
		resp1 <- c(rep(0, delay), resp)[1:l]
		for (i in 1:(n-1)){
			resp <- resp + resp1
			resp1 <- c(rep(0, delay), resp1)[1:l]			
		}
		resp <- resp[1:L]	
	} else {
		if (is.null(ampls)) ampls <- rep(1, length(times))
		l <- length(resp)
		if (L > l) stop ("L is too large")
		times <- round(times / dt, 0)
		if ((max(times)) >= l) stop ("delay is too large")
		resp1 <- resp; resp <- rep(0, l)
		for (i in 1:length(times)){
			resp2 <- ampls[i] * c(rep(0, times[i]), resp1)[1:l]
#			print(ampls[i])			
			resp <- resp + resp2
#			if (i==1) plot(resp, t="l", ylim=c(0,20)) else lines(resp, col=2)
		}
		resp <- resp[1:L]	
	}
	resp
}

comb.lin.mat <- function(resp, delay, dt=0.1){
	if (is.vector(resp)) resp <- matrix(resp, ncol=1)
	L <- nrow(resp)
	resp <- resp - resp[1,1]
	if (is.null(delay)) stop("delay must be specified")
	n <- ncol(resp)
	delay <- delay/dt
	print(c(delay, n))
	if (delay * (n-1) >= L) stop ("delay is too large")
	if (n == 1) return(resp[,1])
	resp1 <- resp[,1]
	for (i in 2:n){
		resp1 <- resp1 + c(rep(0, (i-1) * delay), resp[,i])[1:L]
	}
	resp1
}



# pred.lin <- comb.lin(single, times=times, ampls=seq(1, 20), dt=0.5)
# plot(pred.lin)
# abline(v=(times+12.1)*2)
# t <- seq(0,100, length=1000)
# r <- exp(-t/10); r[1:10] <- 0
# plot(t, r, t="l")

# tt <- runif(20, 0, 90)
# r2 <- comb.lin(r, 5, 5)
# plot(t, r2, t="l")
# r2 <- comb.lin(r, 5, 5, L=500)
# plot(t[1:500], r2, t="l")

# r2 <- comb.lin(r, times=tt)
# plot(t, r2, t="l")
# abline(v=tt+1)
# abline(v=tt[1]+1, col=2)

# r2 <- comb.lin(r, times=tt, ampls=runif(20))
# lines(t, r2, t="l", col=3)
