## Plot the inputs of the simulated L2/3 cells

##########################################

Tmax <- 24
Erate <- 5
Irate <- 40
types <- c('balanced', 'random_NR')
i.type <- 2
type <- types[i.type] 

rseed <- 1
rates <- matrix(NA, 16, 24000)

for (irep in 1:16){ # takes about 2 min on a Macbook with 2.2GHz i7 processor
	#rseed <- 2
	ename <- paste(type, '/Espikes_d', Tmax, '_Ne1920_Re', Erate, '_rseed', rseed, '_rep', irep-1, '.dat', sep='')
	stim <- read.table(ename)
	iname <- paste(type, '/Ispikes_d', Tmax, '_Ni192_Ri', Irate, '_rseed', rseed, '_rep', irep-1, '.dat', sep='')
	istim <- read.table(iname)
	
	dt <- 1/1000		
	spikes <- stim
	# spikes <- istim

	poprate <- rep(0, Tmax/dt)
	spt <- spikes
	for (i in 1:nrow(spt)) {
		poprate[spt[i,2]] <- poprate[spt[i,2]] + 1
	}
		
	sdfilt <- 0.1 # s
	filt <- dnorm(seq(-4, 4, length=sdfilt * 2 * 1000 * 8 + 1))
	filt <- filt / sum(filt)
	filtrate.theta <- filter(poprate, filt, circular = T)
	rates[irep,] <- filtrate.theta
	# matplot(cbind(poprate, filtrate.gamma, filtrate.theta), t='l', col=c(1,2,3), lwd=c(1,2,3))
	# matplot(cbind(filtrate.gamma, filtrate.theta), t='l', col=c(2,3), lwd=c(1,2,3))
	cat(irep, ' ')
}



ii <- seq(1, 24000, by=100)
pdf(file=paste('input_spcounts_', type, '.pdf', sep=''), 8, 3, useDingbats=F)
par(mar=c(4,4,1,1))
matplot(t(rates[,ii]), t='l', lty=1, col=rainbow(16), axes=F, xlab='time (s)', ylab='spike counts / ms', ylim=c(7.8, 8.65))
lines(colMeans(rates[,ii]), t='l', lwd=2)
axis(1, 0:8 * 30, 0:8 * 3)
axis(2, las=2)
dev.off()

ie.cols <- c(rgb(1, 204/255, 153/255), rgb(1, 0, 0), rgb(153/255, 204/255, 1))
outname <- paste('inputspikes_clust_', type, '_d', Tmax, '_Re', Erate, '_Ri', Irate, '_rseed', rseed, '_rep', irep, '.png', sep='')
png(filename=outname, 1600, 600, pointsize=24)
par(mar=c(4,4,2,1))
layout(matrix(c(1,2), 1), c(3,1))
cols <- rep(ie.cols[1], nrow(stim))
cols[(stim[,1] > 1560) & (stim[,1] < 1800)] <- ie.cols[2]
plot(stim[,2], stim[,1], pch=16, cex=0.2, col=cols, axes=F, ylim=c(0, 2200), xlab="time (s)", ylab="cells"); axis(2, las=2)
axis(1, 0:8 * 3000, 0:8 * 3)
points(istim[,2], 1920 + istim[,1], pch=16, cex=0.2, col=ie.cols[3])

t0 <- 6000; t1 <- 7500
ie <- which((stim[,2] > t0) & (stim[,2] < t1))
se <- stim[ie,]
ii <- which((istim[,2] > t0) & (istim[,2] < t1))
si <- istim[ii,]
cols <- rep(ie.cols[1], nrow(se))
cols[(se[,1] > 1560) & (se[,1] < 1800)] <- ie.cols[2]
par(mar=c(4,0,2,1))
plot(se[,2], se[,1], pch=16, cex=0.2, col=cols, axes=F, ylim=c(0, 2200), xlab="time (s)", ylab="")
axis(1, c(t0, t1), c(t0, t1)/1000)
points(si[,2], 1920 + si[,1], pch=16, cex=0.2, col=ie.cols[3])

dev.off()
			