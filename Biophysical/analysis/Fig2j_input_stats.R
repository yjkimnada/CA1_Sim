####################################################################################
## evaluate the variability of the inputs to the biophysical model
####################################################################################

localdir <- getwd()
library(viridis)
library(gplots)

##################################################################
## analyse low level data - skip this part, if you don't have the stimulus files generated
##################################################################
## stimuli should be generated via https://bitbucket.org/bbu20/popact/
## stimulus files can be rather large (spikes of 2000 cells for 10 s), so they are not provided 

datadir <- 'CA1stims/' # change to the local datapath, where the stimulus files can be found
setwd(datadir)

Tmax <- 10
Erate <- 1
Irate <- 7.5
types <- c('balanced', 'random_NR')
sds <- array(0, dim=c(2, 6, 10), dimnames=list(types, c('SE-E', 'sd of inrate-E', 'SE-I', 'sd of inrate-I', 'SE-EI', 'sd of inrate-EI'), paste('rep', 1:10)))
i.type <- 1
rseed <- 5
irep <- 9

for (i.type in 1:2){ # balanced or random
	type <- types[i.type] 
	setwd(type)
	for (rseed in 1:10){ # all random seeds
		#irep <- 1
		erates <- matrix(NA, 16, 10000)
		irates <- matrix(NA, 16, 10000)
		rates <- matrix(NA, 16, 10000)
		for (irep in 1:16){ # all repetitions
			#rseed <- 2
			## excitatory population activity - cell# and time (ms)
			ename <- paste('Espikes_d', Tmax, '_Ne2000_Re', Erate, '_rseed', rseed, '_rep', irep-1, '.dat', sep='')
			stim <- read.table(ename)

			## excitatory population activity - cell# and time (ms)
			iname <- paste('Ispikes_d', Tmax, '_Ni200_Ri', Irate, '_rseed', rseed, '_rep', irep-1, '.dat', sep='')
			istim <- read.table(iname)
						
			dt <- 1/1000		
			spikes <- stim
			ispikes <- istim

			poprate <- rep(0, Tmax/dt)
			ipoprate <- rep(0, Tmax/dt)

			for (i in 1:nrow(spikes)) {
				poprate[spikes[i,2]] <- poprate[spikes[i,2]] + 1
			}
			poprate <- poprate[1:(Tmax/dt)]
			for (i in 1:nrow(ispikes)) {
				ipoprate[ispikes[i,2]] <- ipoprate[ispikes[i,2]] + 1
			}
			ipoprate <- ipoprate[1:(Tmax/dt)]
			
			sdfilt <- 0.1 # s
			filt <- dnorm(seq(-4, 4, length=sdfilt * 2 * 1000 * 8 + 1))
			filt <- filt / sum(filt)
			filtErate.theta <- filter(poprate, filt, circular = T)
			filtIrate.theta <- filter(ipoprate, filt, circular = T)
			filtrate.theta <- filter(poprate - 2 * ipoprate, filt, circular = T)

			erates[irep,] <- filtErate.theta
			irates[irep,] <- filtIrate.theta
			rates[irep,] <- filtrate.theta

			cat(irep, ' ')
		}
		mrate <- colMeans(erates)
		sds[i.type, 2, rseed] <- sd(mrate)
		sds[i.type, 1, rseed] <- mean(apply(erates, 2, sd)) / sqrt(16)

		mrate <- colMeans(irates)
		sds[i.type, 4, rseed] <- sd(mrate)
		sds[i.type, 3, rseed] <- mean(apply(irates, 2, sd)) / sqrt(16)

		mrate <- colMeans(rates)
		sds[i.type, 6, rseed] <- sd(mrate)
		sds[i.type, 5, rseed] <- mean(apply(rates, 2, sd)) / sqrt(16)

		cat('\n seed:', rseed, ' \n')
	}
	cat('\n type:', type, ' \n')
	setwd('../')
}

setwd(localdir)
save(sds, file='./datasets/sds_inspikes.Rdata')

##################################################################
## generate Figure from the high level data
##################################################################
setwd(localdir)
load('./datasets/sds_inspikes.Rdata')

Pbal1 <- t.test(sds[1,1,] - sds[1,2,])$p.value # balanced Excitatory input differs from baseline?
Pbal2 <- t.test(sds[1,3,] - sds[1,4,])$p.value # balanced Inhibitory input differs from baseline?
Pbal3 <- t.test(sds[1,5,] - sds[1,6,])$p.value # balanced total input differs from baseline?


##############################################
cols <- c(rgb(255, 204, 170, max=255), rgb(140, 208, 255, max=255), rgb(200, 170, 120, max=255))

# pdf(file='Fig2j_input_stats.pdf', 3, 3, useDingbats=F)
par(mar=c(5,5,2,1))

mm <- rowMeans(sds[c(2,1),2,]^2)
sems <- apply(sds[c(2,1),2,]^2, 1, sd) / sqrt(10)
plotCI(1:2+0.15, mm, sems, ylim=c(0, 0.02), xlim=c(0, 3), main='excitatory', axes=F, pt.bg=cols[1], pch=22, gap=0, ylab='variance (kHz^2)', xlab='')
matplot(1:2-0.15, sds[c(2,1),2,]^2, pch=16, cex=0.5, col=grey(0.75), add=T)
segments(x0=1:2-0.2,x1=1:2+0.2, rowMeans(sds[c(2,1),1,]^2), col=1, lwd=1)
segments(x0=1:2-0.2,x1=1:2+0.2, rowMeans(sds[c(2,1),1,]^2)*16, col=1, lwd=1)
axis(1, 1:2, c('random', 'uniform'), tick=F, las=2)
axis(2, las=2)
# dev.off()

write.table(format(round(sds[c(2,1),2,]^2, 6), scientific=F), file='Fig2j_data.txt')

