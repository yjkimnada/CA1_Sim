## Plot the inputs of the simulated HPC place cells

##########################################

types <- c('balanced', 'random_NR')
type <- types[2] 
setwd(type)

irep <- 1
rseed <- 1
Tmax <- 10
Erate <- 1
Irate <- 7.5

ename <- paste('Espikes_d', Tmax, '_Ne2000_Re', Erate, '_rseed', rseed, '_rep', irep, '.dat', sep='')
stim <- read.table(ename)

iname <- paste('Ispikes_d', Tmax, '_Ni200_Ri', Irate, '_rseed', rseed, '_rep', irep, '.dat', sep='')
istim <- read.table(iname)

setwd('../')


ie.cols <- c(rgb(239/255, 147/255, 0), rgb(20/255, 204/255, 0), rgb(0, 154/255, 204/255))

outname <- paste('inputspikes_clust_', type, '_d', Tmax, '_Re', Erate, '_Ri', Irate, '_rseed', rseed, '_rep', irep, '.png', sep='')
png(filename=outname, 1200, 600, pointsize=24)
cols <- rep(ie.cols[1], nrow(stim))
cols[(stim[,1] > 880) & (stim[,1] < 1120)] <- ie.cols[2]
plot(stim[,2], stim[,1], pch=16, cex=0.2, col=cols, axes=F, ylim=c(0, 2200), xlab="time (ms)", ylab="cells"); axis(1); axis(2, las=2)
points(istim[,2], 2000 + istim[,1], pch=16, cex=0.2, col=ie.cols[3])
dev.off()
			
