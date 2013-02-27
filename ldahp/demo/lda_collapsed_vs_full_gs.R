library(plotrix);
library(ldahp);

## Initialize variables
K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 20000 # the maximum number of Gibbs iterations
burn.in        <- 8000
spacing        <- 1
lambda.hat     <- 80
gen.alpha.v    <- c(3, 2) # symmetric Dirichlet
eta            <- 1

set.seed(1983)

## Generates documents with a given beta.m
beta.m              <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, 1:12]     <- array(10, c(1, 12))
beta.m[2, 8:20]     <- array(10, c(1, 13))
ds                  <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);


## The full Gibbs sampling
ptm                <- proc.time();
fg.mdl             <- lda_full_c2(K, V, ds, gen.alpha.v, eta, max.iter, burn.in, spacing);
ptm                <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling
ptm               <- proc.time();
cg.mdl            <- lda_collapsed_gibbs_c(K, V, ds, gen.alpha.v, eta, max.iter, burn.in, spacing);
ptm               <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


par(mfrow = c(2, 2))

# log marginal posterior plots 

x.axis <- 1:dim(fg.mdl$lmp)[1]
plot(x.axis, fg.mdl$lmp[,1], type="l", col="blue", ylab="log marginal posterior", xlab="Gibbs iterations", main="Full GS", lwd=0.4, ylim=c(60800, 62200)) 
plot(x.axis, cg.mdl$lmp[,1], type="l", col="black", ylab="log marginal posterior", xlab="Gibbs iterations", main="Collapsed GS", lwd=0.4, ylim=c(60800, 62200)) 

# acf on the log marginal posterior 

acf(fg.mdl$lmp[,1], lag.max=100, main=expression(paste("Full GS: log marginal posterior(acf)")))
acf(cg.mdl$lmp[,1], lag.max=100, main=expression(paste("Collapsed GS: log marginal posterior(acf)")))



calc_mean_sd <- function(start.idx, end.idx, lmp){
  slmp <- lmp[start.idx:end.idx]
  list(mean=mean(slmp), sd=sd(slmp))
}


calc_mean_sd(8000, 9000, fg.mdl$lmp[,1])
calc_mean_sd(9000, 10000, fg.mdl$lmp[,1])
calc_mean_sd(10000, 11000, fg.mdl$lmp[,1])
calc_mean_sd(11000, 12000, fg.mdl$lmp[,1])


calc_mean_sd(8000, 9000, cg.mdl$lmp[,1])
calc_mean_sd(9000, 10000, cg.mdl$lmp[,1])
calc_mean_sd(10000, 11000, cg.mdl$lmp[,1])
calc_mean_sd(11000, 12000, cg.mdl$lmp[,1])



# Autocorrelation plots 

# acf computed on the theta counts 

par(mfrow = c(3,2))
acf(fg.mdl$theta[1,1,], lag.max=100, main=expression(paste("Full GS: ", theta)[1][1]))
acf(cg.mdl$theta[1,1,], lag.max=100, main=expression(paste("Collapsed GS: ", theta)[1][1]))

acf(fg.mdl$theta[1,8,], lag.max=100, main=expression(paste("Full GS: ", theta)[8][1]))
acf(cg.mdl$theta[1,8,], lag.max=100, main=expression(paste("Collapsed GS: ", theta)[8][1]))

acf(fg.mdl$theta[1,4,], lag.max=100, main=expression(paste("Full GS: ", theta)[4][1]))
acf(cg.mdl$theta[1,4,], lag.max=100, main=expression(paste("Collapsed GS: ", theta)[4][1]))

# acf computed on the beta counts 

par(mfrow = c(3,2))
acf(fg.mdl$beta[1,1,], lag.max=100, main=expression(paste("Full GS: ", beta)[1][1]))
acf(cg.mdl$beta[1,1,], lag.max=100, main=expression(paste("Collapsed GS: ", beta)[1][1]))

acf(fg.mdl$beta[1,8,], lag.max=100, main=expression(paste("Full GS: ", beta)[8][1]))
acf(cg.mdl$beta[1,8,], lag.max=100, main=expression(paste("Collapsed GS: ", beta)[8][1]))

acf(fg.mdl$beta[1,5,], lag.max=100, main=expression(paste("Full GS: ", beta)[5][1]))
acf(cg.mdl$beta[1,5,], lag.max=100, main=expression(paste("Collapsed GS: ", beta)[5][1]))

# Beta plots for the two types of Gibbs samplers 

num.samples <- dim(fg.mdl$beta)[3]
fg.beta <- matrix(0, nrow=K, ncol=V)
cg.beta <- matrix(0, nrow=K, ncol=V)
for (i in (num.samples-100):num.samples){
  fg.beta <- fg.beta + fg.mdl$beta[,,i] / 100;
  cg.beta <- cg.beta + cg.mdl$beta[,,i] / 100;
}

par(mfrow = c(3,1))
color2D.matplot(beta.m, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Data generating ", beta, " matrix")));
color2D.matplot(fg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Full GS output ", beta, " matrix")));
color2D.matplot(cg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Collapsed GS output ", beta, " matrix")));


## Saves every object into a file
rdata_file     <- "ldac_full_collapsed_gibbs.RData"
save.image(rdata_file)

