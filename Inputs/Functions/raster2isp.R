
library("Matrix")

## single cell raster to spike time
raster2isp <- function(raster){
	isp <- which(raster>0)
	if (length(isp) > 0){
		spt <- rep(isp[1], raster[isp[1]])
		isp <- isp[-1]
		if (length(isp) > 0){
			for (i in isp)	spt <- c(spt, rep(i, raster[i]))
		}
	} else spt<- NULL
	spt
}


## raster matrix to spike times
raster2sp <- function(raster){
	if (sum(raster) > 0){
		spt <- c(0,0)
		for (i in 1:nrow(raster)){
			if (sum(raster[i,]) > 0){
				isp <- raster2isp(raster[i,])
				spt <- rbind(spt, cbind(rep(i, length(isp)), isp))
			}
		}
		spt <- spt[-1,]
		i.spt <- sort(spt[,2], index=T)$ix
		spt <- spt[i.spt,]
	} else spt <- NULL
	spt
}


### plot the average ISI histogram of cells in the raster
listRaster2hist <- function(rast, max.ISI=999, dISI=5, add=F, ...){
	N.cells <- nrow(rast[[1]])
	br <- c(seq(0, max.ISI, by=dISI), Inf)
	br.plot <- seq(dISI/2, length=length(br)-1, by=dISI)
	nn <- length(br)
	ISIs <- matrix(0, N.cells, nn-1)
	
	for (i.trial in 1:length(rast)){
		spt <- raster2sp(rast[[i.trial]])
		for (i.cell in 1:N.cells){
			sp <- spt[spt[,1]==i.cell, 2]
			ISIs[i.cell,] <- ISIs[i.cell,] + hist(diff(sp), br=br, plot=F)$counts
		}	
	}
	if (add)	lines(br.plot[1:(nn-2)], colSums(ISIs)[1:(nn-2)], t="s", lty=1, ...) else {
		plot(br.plot[1:(nn-2)], colSums(ISIs)[1:(nn-2)], t="s", lty=1, xlab="time (timebins)", ylab="ISI count", axes=F, ...); axis(1); axis(2, las=2)
	}
}

## convert spike times 
spt2raster <- function(spt, dt, Tmin=NULL, Tmax=NULL){
	if (!(colnames(spt)[1] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	if (!(colnames(spt)[2] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	
	N.cells <- length(unique(spt[,"cell"]))
	ii.cells <- unique(spt[,"cell"])

	if (is.null(Tmin)) Tmin <- min(spt[,"time"])
	if (is.null(Tmax)) Tmax <- max(spt[,"time"])
	spt[,"time"] <- spt[,"time"] - Tmin
	spt[spt[,"time"]==0,"time"] <- dt/2

	Tmax <- Tmax - Tmin
	
	Lmax <- floor(Tmax / dt) + 1
	NN <- Lmax * N.cells
	N.spikes <- nrow(spt)
	sparseness <- N.spikes / NN
	if (sparseness < 1/50)	rast <- Matrix(0, Lmax, N.cells)	else rast <- matrix(0, Lmax, N.cells)
	for (i in 1:N.spikes){
		i.sp <- ceiling(spt[i,"time"] / dt)
		i.cell <- which( ii.cells == spt[i,"cell"])
		rast[i.sp, i.cell] <- rast[i.sp, i.cell] + 1
	}
	rast
}

