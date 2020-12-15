####################################################################################
## A script to analyse the effect of different factors on the membrane potential
####################################################################################

library(viridis)
library(colormap)
library(gplots)

source('./functions/StatsEval_functions.R', chdir=T)

stim_types <- c('random_NR', 'balanced')
act_types <- c('Dactive', 'Dactive')
spinetype <- 'spines'
typenames <- c('act_rand', 'act_bal', 'pass_rand', 'pass_bal')
n.sim <- length(stim_types)

ntrial <- 16
Tmax <- 10
gN <- '0.8'

graphics <- F

## to obtain this dataset, run the script 'functions/Collect_VmStats.R' on the dataset with locally balanced connectivity (inputs: both random and uniform)
load('datasets/response_variability_ttt_globalRand_localReg_spines.Rdata')

sds_L <- stats_data$sds
meanresps_L <- stats_data$meanresps
ttt_variance_L <- stats_data$ttt_variance
ttt_variance_peak_L <- stats_data$ttt_variance_peak

load('datasets/response_variability_ttt_globalRand_localRand_spines.Rdata')

sds_R <- stats_data$sds
meanresps_R <- stats_data$meanresps
ttt_variance_R <- stats_data$ttt_variance
ttt_variance_peak_R <- stats_data$ttt_variance_peak

load('datasets/response_variability_ttt_globalReg_localReg_spines.Rdata')

sds_G <- stats_data$sds
meanresps_G <- stats_data$meanresps
ttt_variance_G <- stats_data$ttt_variance
ttt_variance_peak_G <- stats_data$ttt_variance_peak

load('datasets/response_variability_ttt_globalReg_localRand_spines.Rdata')

sds_Glr <- stats_data$sds
meanresps_Glr <- stats_data$meanresps
ttt_variance_Glr <- stats_data$ttt_variance
ttt_variance_peak_Glr <- stats_data$ttt_variance_peak


############################################################
## all active simulations
require(abind)
sds <- abind(sds_R, sds_L, sds_Glr, sds_G, along=2)
sds_act <- sds[,c(1,3,5,7,2,4,6,8),]

dimnames_sd <- dimnames(sds_act)
dimnames_sd[[2]] <- c('rand_rand', 'rand_loc', 'rand_glob', 'rand_full', 'unif_rand', 'unif_loc', 'unif_glob', 'unif_full')
dimnames(sds_act) <- dimnames_sd

ttt_variance <- abind(ttt_variance_R, ttt_variance_L, ttt_variance_Glr, ttt_variance_G, along=1)
ttt_variance_act <- ttt_variance[c(1,3,5,7,2,4,6,8),,]

ttt_variance_peak <- abind(ttt_variance_peak_R, ttt_variance_peak_L, ttt_variance_peak_Glr, ttt_variance_peak_G, along=2)
ttt_variance_peak_act <- ttt_variance_peak[,c(1,3,5,7,2,4,6,8),]

typenames <- paste(rep(c('iR', 'iU'), e=4), c('sR', 'sL', 'sG', 'sF'))

col_iR_sR <- rgb(200, 183, 238, max=255)
col_iR_sL <- rgb(226, 217, 246, max=255)
col_iR_sG <- rgb(255, 161, 185, max=255)
col_iR_sF <- rgb(255, 208, 220, max=255)
col_iU_sR <- rgb(140, 208, 254, max=255)
col_iU_sL <- rgb(210, 236, 255, max=255)
col_iU_sG <- rgb(181, 181, 181, max=255)
col_iU_sF <- rgb(231, 231, 231, max=255)

cols <- c(col_iR_sR, col_iR_sL, col_iR_sG, col_iR_sF, col_iU_sR, col_iU_sL, col_iU_sG, col_iU_sF)
pchs <- c(22, 25, 22, 25, 22, 25, 22, 25)

#############################################
## Fig 2k
#############################################

## three different way to plot the same thing: variability of row membrane potential, SD and variance of the filtered Vm:
fname <- 'fig2k_response_variability.pdf'
plot_data_SEM(filename=NULL, sds_act, ttt_variance_act, ttt_variance_peak_act, typenames, col=cols, pchs=pchs, plotvars=T)
# plot_data_SEM(filename=fname, sds_act, ttt_variance_act, ttt_variance_peak_act, typenames, col=cols, pchs=pchs, plotvars=T)

write.table(round(sds_act[1,,],4), file='Fig2k_data.txt')

#############################################
### statistical significance
wilcox.test(sds_act[1,1,]^2, sds_act[1,2,]^2, paired=TRUE, conf.int=T) # 0.38
wilcox.test(sds_act[1,1,]^2, sds_act[1,4,]^2, paired=TRUE, conf.int=T) # 0.03
wilcox.test(sds_act[1,3,]^2, sds_act[1,4,]^2, paired=TRUE, conf.int=T) # 0.5

wilcox.test(sds_act[1,1,]^2, sds_act[1,5,]^2, paired=TRUE, conf.int=T) # 0.002
wilcox.test(sds_act[1,5,]^2, sds_act[1,6,]^2, paired=TRUE, conf.int=T) # 0.77
wilcox.test(sds_act[1,5,]^2, sds_act[1,8,]^2, paired=TRUE, conf.int=T) # 0.002
wilcox.test(sds_act[1,7,]^2, sds_act[1,8,]^2, paired=TRUE, conf.int=T) # 0.23
wilcox.test(sds_act[1,8,]^2, apply(ttt_variance_act[8,,], 1, mean) / 16, paired=TRUE, conf.int=T) # 0.77

#############################################
## Fig 2j
#############################################

###############################################################
i_decomp <- c(1,7, 2,8,1) # iR_sR, iR_sF, iU_sR, iU_sU
sds_decomp <- sds[,i_decomp,]
ttt_variance_decomp <- ttt_variance[i_decomp,,]
ttt_variance_peak_decomp <- ttt_variance_peak[, i_decomp,]
typenames <- c('random', 'uniform',  'regular', 'uni-reg', 'expected')

ttt_Var_decomp <- apply(ttt_variance_decomp/16, c(1,2), mean)[3,]
sds_decomp[1,5,] <- sqrt(sds_decomp[1,2,]^2 + sds_decomp[1,3,]^2 - ttt_Var_decomp)
ttt_Var_peak_decomp <- apply(ttt_variance_peak_decomp, c(2,3), mean)[3,]
sds_decomp[3,5,] <- sds_decomp[3,2,] + sds_decomp[3,3,] - ttt_Var_peak_decomp
ttt_variance_decomp[5,,] <- NA
ttt_variance_peak_decomp[,5,] <- NA


typenames <- c('dendritic', 'input', 'dend & inp', 'dend + inp', 'clustering-rI',  'clustering-uI')

ttt_Var_dend <- apply(ttt_variance_decomp/16, c(1,2), mean)[3,]
tuning_var.dend <- sds_decomp[1,3,]^2 - ttt_Var_dend

ttt_Var_inp <- apply(ttt_variance_decomp/16, c(1,2), mean)[2,]
tuning_var.inp <- sds_decomp[1,2,]^2 - ttt_Var_inp

ttt_Var_di <- apply(ttt_variance_decomp/16, c(1,2), mean)[1,]
tuning_var.di <- sds_decomp[1,1,]^2 - ttt_Var_di


tuning_vars <- cbind(tuning_var.dend, tuning_var.inp, tuning_var.di)
mm <- colMeans(tuning_vars)
ss <- apply(tuning_vars, 2, sd) / sqrt(10)

mm_exp <- mm[1] + mm[2]
ss_exp <- sqrt(var(tuning_vars[,1]) + var(tuning_vars[,2])) / sqrt(10)

## we estimate the effect of clustering by substracting the locally balanced from the random - in both the random input case (not shown in the paper)
mm_sc1 <- mean(sds_act[1,1,] - sds_act[1,2,])
ss_sc1 <- sqrt(var(sds_act[1,2,]) + var(sds_act[1,1,])) / sqrt(10)

## and also in the balanced input case
mm_sc2 <- mean(sds_act[1,5,] - sds_act[1,6,])
ss_sc2 <- sqrt(var(sds_act[1,5,]) + var(sds_act[1,6,])) / sqrt(10)

mm <- c(mm, mm_exp, mm_sc1, mm_sc2)
ss <- c(ss, ss_exp, ss_sc1, ss_sc2)


# pdf(file='fig2j_factor_variances.pdf', 6, 3, useD=F)
par(mar=c(6,4,2,1))

plotCI(1:6+0.15, mm, ss, gap=0, pch=22, lwd=1, pt.bg=cols, cex=2, main='factor variances', xlab='', ylab='voltage2 (mV^2)', axes=F, xlim=c(0,7), ylim=c(-0.025, 0.1))
matplot(1:3-0.15, t(tuning_vars[,1:3]), col=grey(0.74), pch=16, cex=0.7, add=T)	
abline(h=0)
axis(2, las=2)
axis(1, 1:6, typenames, las=2, tick=F)
# dev.off()
