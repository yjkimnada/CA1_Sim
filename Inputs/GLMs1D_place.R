## Script to generate hippocampal population activity during theta

source('./Functions/gen_cells.R', chdir=T) # load the required functions - essentially simulating GLMs in R

Ensyn <- 2000; Insyn=200; Ntrial <- 1000; Erate <- 0.5; Irate=7.4; dt <- 0.001; Erate.bg <- 0.1 # these are the original parameters
#Ensyn <- 2000; Insyn=200; Ntrial <- 651; Erate <- 1; Irate=7.5; dt <- 0.001; Erate.bg <- 1
types <- c("balanced", "random_NR")
rseed <- 1 # set rseed to the random seed 1-10
itype <- 2 # set itype to 1 (balanced) or 2 (random)

graphics <- F

### we change a few parameters to speed up the demo
### comment this out if you want to replicate the original inputs
#Ensyn <- 40 
#Insyn <- 10
#Ntrial <- 4


# loading w.template for phase precession - fitted to Skaggs et al., 1996
infile <- paste('./Functions/wTemplateSkaggs_L50nx10nphi4.RData', sep='')
load(file=infile)

type <- types[itype]

sample.template=F; randfields=F; randrates=F; fields.middle=F; same.runs=T
if (type == "random_N") randfields <- T
if (type == "random_NR") {
	randfields <- T
	randrates <- T
}


sp.E <- gen.placecells(N.cells=Ensyn, L =200, Ntrial= Ntrial, dt=dt, N.x.basis=40, N.phi.basis=4, rate.ave=Erate, rate.bg=Erate.bg, seed=rseed, mu.v=0.2, vx=NULL, phi=NULL, graphics=graphics, plot.phase.precess=F, adapt=T, sample.template=sample.template, randfields=randfields, randrates=randrates, fields.middle=fields.middle, var.trials=1/100000, same.runs=same.runs, w.template=w.template)
		
cat("inhibitory simulations started... \n")
	
sp.I <- gen.thetacells(N.cells=Insyn, L =100, Ntrial= Ntrial, dt=dt, N.phi.basis=4, rate.ave=Irate, seed= rseed + 100, vx=list(x=sp.E$dists, v=sp.E$speed.1D), phi=sp.E$lfp[[1]], graphics=graphics, adapt=T, same.cells=T)
		
outdir <- paste("./Data/place/", type, sep="")
dir.create(file.path(outdir), showWarnings = TRUE)
duration <- round(length(sp.E$dists[[1]]) * dt)
		
for (rep in 1:Ntrial){

	Espikes <- sp.E$spt[[rep]][,c(2,1)] # first column: cell, second column: time in ms
	Espikes[,1] <- Espikes[,1] - 1
	Espikes[,2] <- (Espikes[,2] - sp.E$t[rep, 1]) * 1000
	fname  <- paste (outdir, '/Espikes_d', duration, '_Ne', Ensyn, '_Re', Erate, '_rseed', rseed, '_rep', rep-1, '.dat', sep="")
	write.table(Espikes, file=fname, row.names=FALSE, col.names=FALSE)

	Ispikes <- sp.I$spt[[rep]][,c(2,1)] # first column: cell, second column: time in ms
	Ispikes[,1] <- Ispikes[,1] - 1
	Ispikes[,2] <- round((Ispikes[,2] - sp.I$t[rep, 1]) * 1000)
	fname  <- paste (outdir, '/Ispikes_d', duration, '_Ni', Insyn, '_Ri', Irate, '_rseed', rseed, '_rep', rep-1, '.dat', sep="")
	write.table(Ispikes, file=fname, row.names=FALSE, col.names=FALSE)
 	
 }

