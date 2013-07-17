##############################################################################################
## This script is used to generate comparison plots of the LDA collapsed 
## and full Gibbs samplers' log marginal posterior on z, \theta, and \beta 
## 
## Created on: July 13, 2013 
##############################################################################################

# Init file names 

rdata.file <- "/home/clintpg/results/fg_cg_ae33.RData"
log.marginal.posterior.plot <- "/home/clintpg/results/fg_cg_ae33_lmp.eps"
theta.acf.plot <- "/home/clintpg/results/fg_cg_ae33_theta.eps"
beta.acf.plot <- "/home/clintpg/results/fg_cg_ae33_beta.eps"


# Loads the saved data 

load(rdata.file)


# --------------------------------------------------------------------------------------
# Log marginal posterior plots 
# --------------------------------------------------------------------------------------

postscript(file=log.marginal.posterior.plot, title="log marginal posterior plots", horiz=F)  
par(mfrow = c(2, 2))

x.axis <- (burn.in+1):(burn.in+dim(fg.mdl$lmp)[1])
plot(x.axis, fg.mdl$lmp[,1], type="l", col="blue", ylab="log marginal posterior", xlab="GS iterations", main=expression(paste("GS on (", beta, ", ", theta, ", z): log marginal posterior (trace)")), lwd=0.4, ylim=c(59500, 60300), cex = 1.5, cex.lab = 1.5, cex.main = 1.2) 
plot(x.axis, cg.mdl$lmp[,1], type="l", col="black", ylab="log marginal posterior", xlab="GS iterations", main=expression(paste("GS on (z): log marginal posterior (trace)")), lwd=0.4, ylim=c(59500, 60300), cex = 1.5, cex.lab = 1.5, cex.main = 1.2) 

# acf on the log marginal posterior 

acf(fg.mdl$lmp[,1], lag.max=100, main=expression(paste("GS on (", beta, ", ", theta, ", z): log marginal posterior (acf)")), cex = 1.5, cex.lab = 1.5, cex.main = 3)
acf(cg.mdl$lmp[,1], lag.max=100, main=expression(paste("GS on (z): log marginal posterior (acf)")), cex = 1.5, cex.lab = 1.5, cex.main = 3)

dev.off()



# --------------------------------------------------------------------------------------
# Autocorrelation plots 
# --------------------------------------------------------------------------------------

# acf computed on the theta counts 

postscript(file=theta.acf.plot, title="theta elements' trace plots", horiz=F)  
par(mfrow = c(2,2))
acf(fg.mdl$theta[1,1,], lag.max=50, main=expression(paste("GS on (", beta, ", ", theta, ", z): sample of ", theta)[1][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
acf(cg.mdl$theta[1,1,], lag.max=50, main=expression(paste("GS on (z): estimate of ", theta)[1][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)

acf(fg.mdl$theta[1,8,], lag.max=50, main=expression(paste("GS on (", beta, ", ", theta, ", z): sample of ", theta)[8][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
acf(cg.mdl$theta[1,8,], lag.max=50, main=expression(paste("GS on (z): estimate of ", theta)[8][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
dev.off()

# acf computed on the beta counts 

postscript(file=beta.acf.plot, title="beta elements' trace plots", horiz=F)  
par(mfrow = c(2,2))
acf(fg.mdl$beta[1,1,], lag.max=50, main=expression(paste("GS on (", beta, ", ", theta, ", z): sample of ", beta)[1][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
acf(cg.mdl$beta[1,1,], lag.max=50, main=expression(paste("GS on (z): estimate of ", beta)[1][1]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)

acf(fg.mdl$beta[1,4,], lag.max=50, main=expression(paste("GS on (", beta, ", ", theta, ", z): sample of ", beta)[1][4]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
acf(cg.mdl$beta[1,4,], lag.max=50, main=expression(paste("GS on (z): estimate of ", beta)[1][4]), lwd=3, cex = 1.5, cex.lab = 1.5, cex.main = 1.5)
dev.off()





# # --------------------------------------------------------------------------------------
# # Beta plots for the two types of Gibbs samplers 
# # --------------------------------------------------------------------------------------
# 
# num.samples <- dim(fg.mdl$beta)[3]
# fg.beta <- matrix(0, nrow=K, ncol=V)
# cg.beta <- matrix(0, nrow=K, ncol=V)
# for (i in (num.samples-100):num.samples){
#   fg.beta <- fg.beta + fg.mdl$beta[,,i] / 100;
#   cg.beta <- cg.beta + cg.mdl$beta[,,i] / 100;
# }
# 
# library(plotrix);
# par(mfrow = c(3,1))
# color2D.matplot(beta.m, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Synthetic data generating ", beta, " matrix")));
# color2D.matplot(fg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Full GS: estimated ", beta, " matrix")));
# color2D.matplot(cg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Collapsed GS: estimated ", beta, " matrix")));
# 
# 
# 
# 
# 
# calc_mean_sd <- function(start.idx, end.idx, lmp){
#   slmp <- lmp[start.idx:end.idx]
#   list(mean=mean(slmp), sd=sd(slmp))
# }
# 
# 
# calc_mean_sd(2000, 3000, fg.mdl$lmp[,1])
# calc_mean_sd(3000, 4000, fg.mdl$lmp[,1])
# calc_mean_sd(4000, 5000, fg.mdl$lmp[,1])
# calc_mean_sd(5000, 6000, fg.mdl$lmp[,1])
# 
# calc_mean_sd(2000, 3000, cg.mdl$lmp[,1])
# calc_mean_sd(3000, 4000, cg.mdl$lmp[,1])
# calc_mean_sd(4000, 5000, cg.mdl$lmp[,1])
# calc_mean_sd(5000, 6000, cg.mdl$lmp[,1])
# 
