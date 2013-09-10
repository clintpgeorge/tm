##############################################################################################
## This script is used to generate comparison plots of the LDA collapsed 
## and full Gibbs samplers' log marginal posterior on z, \theta, and \beta 
## 
## Created on: July 13, 2013 
##############################################################################################

## Loads packages 

library(ldahp); 

## Initialize variables

set.seed(1983)

K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 11000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 1
lambda.hat     <- 80
gen.eta        <- 3
gen.alpha      <- 3
gen.alpha.v    <- array(gen.alpha, c(1, K));           # symmetric Dirichlet
gen.eta.v      <- array(gen.eta, c(1, V));             # symmetric Dirichlet
rdata.file     <- "/home/clintpg/results/fg_cg_ea33.RData"
store.Dir      <- 1                                    # store the \theta and \beta Dirichlet samples ? 

## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);



## The full Gibbs sampling

ptm                <- proc.time();
fg.mdl             <- lda_full_c2(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm                <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm               <- proc.time();
cg.mdl            <- lda_collapsed_gibbs_c(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm               <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Saves every object into a file

save.image(rdata.file)


##############################################################################################

# Init file names 

rdata.file <- "/home/clintpg/results/fg_cg_ea33.RData"
log.marginal.posterior.plot <- "/home/clintpg/results/fg_cg_ea33_lmp.eps"
log.posterior.plot <- "/home/clintpg/results/fg_cg_ea33_lp.eps"
theta.acf.plot <- "/home/clintpg/results/fg_cg_ea33_theta.eps"
beta.acf.plot <- "/home/clintpg/results/fg_cg_ea33_beta.eps"


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
# Log posterior plots 
# --------------------------------------------------------------------------------------

postscript(file=log.posterior.plot, title="Log posterior plots", horiz=F)  
par(mfrow = c(2, 2))

x.axis <- (burn.in+1):(burn.in+dim(fg.mdl$lmp)[1])
plot(x.axis, fg.mdl$lmp[,1], type="l", col="blue", ylab="Log posterior", xlab="iteration", main=expression(paste("GS on (", beta, ", ", theta, ", z): log marginal posterior (trace)")), lwd=0.4, ylim=c(59500, 60300), cex = 1.5, cex.lab = 1.5, cex.main = 1.2) 
plot(x.axis, cg.mdl$lmp[,1], type="l", col="black", ylab="Log posterior", xlab="iteration", main=expression(paste("GS on (z): log marginal posterior (trace)")), lwd=0.4, ylim=c(59500, 60300), cex = 1.5, cex.lab = 1.5, cex.main = 1.2) 

# acf on the log marginal posterior 

acf(fg.mdl$lmp[,1], lag.max=100, main=expression(paste("GS on (", beta, ", ", theta, ", z): log posterior (acf)")), cex = 1.5, cex.lab = 1.5, cex.main = 3)
acf(cg.mdl$lmp[,1], lag.max=100, main=expression(paste("GS on (z): log posterior (acf)")), cex = 1.5, cex.lab = 1.5, cex.main = 3)

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
