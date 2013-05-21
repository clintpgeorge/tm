library(plotrix);
library(ldahp);

## Initialize variables

set.seed(1983)

# ## for three topics 
# K              <- 3 # the number of topics
# D              <- 100 # the total number of documents to be generated
# V              <- 40 # the vocabulary size
# max.iter       <- 10000 # the maximum number of Gibbs iterations
# burn.in        <- 4000
# spacing        <- 1
# lambda.hat     <- 80
# gen.alpha.v    <- c(3, 2, 1)
# eta            <- 1
# 
# beta.m              <- matrix(1e-2, nrow=K, ncol=V)
# beta.m[1, 1:12]     <- 10
# beta.m[2, 8:25]     <- 10
# beta.m[3, 15:40]     <- 10 
# rdata_file     <- "ldac_fg_vs_cg_3t.RData"

## for two topics 
K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 10000 # the maximum number of Gibbs iterations
burn.in        <- 4000
spacing        <- 1
lambda.hat     <- 80
gen.alpha.v    <- c(3, 2)
eta            <- 1

beta.m              <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, 1:12]     <- 10
beta.m[2, 8:20]     <- 10

rdata_file     <- "ldac_fg_vs_cg_2t.RData"


## Generates documents with a given beta.m

ds                  <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);



## The full Gibbs sampling
ptm                <- proc.time();
fg.mdl             <- lda_full_c2(K, V, ds, gen.alpha.v, eta, max.iter, burn.in, spacing, 1);
ptm                <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling
ptm               <- proc.time();
cg.mdl            <- lda_collapsed_gibbs_c(K, V, ds, gen.alpha.v, eta, max.iter, burn.in, spacing, 1);
ptm               <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Saves every object into a file

save.image(rdata_file)



par(mfrow = c(2, 2))

# log marginal posterior plots 
# 2 topics ylim=c(60800, 62200)
# 2 topics ylim=c(53200, 55500)

x.axis <- (burn.in+1):(burn.in+dim(fg.mdl$lmp)[1])
plot(x.axis, fg.mdl$lmp[,1], type="l", col="blue", ylab="log marginal posterior", xlab="GS iterations", main=expression(paste("Full GS: samples ", beta, ", ", theta, ", ", z)), lwd=0.4, ylim=c(60800, 62200)) 
plot(x.axis, cg.mdl$lmp[,1], type="l", col="black", ylab="log marginal posterior", xlab="GS iterations", main=expression(paste("Collapsed GS: samples ", z)), lwd=0.4, ylim=c(60800, 62200)) 

# acf on the log marginal posterior 

acf(fg.mdl$lmp[,1], lag.max=100, main=expression(paste("Full GS: log marginal posterior(acf)")))
acf(cg.mdl$lmp[,1], lag.max=100, main=expression(paste("Collapsed GS: log marginal posterior(acf)")))



calc_mean_sd <- function(start.idx, end.idx, lmp){
  slmp <- lmp[start.idx:end.idx]
  list(mean=mean(slmp), sd=sd(slmp))
}


calc_mean_sd(2000, 3000, fg.mdl$lmp[,1])
calc_mean_sd(3000, 4000, fg.mdl$lmp[,1])
calc_mean_sd(4000, 5000, fg.mdl$lmp[,1])
calc_mean_sd(5000, 6000, fg.mdl$lmp[,1])

calc_mean_sd(2000, 3000, cg.mdl$lmp[,1])
calc_mean_sd(3000, 4000, cg.mdl$lmp[,1])
calc_mean_sd(4000, 5000, cg.mdl$lmp[,1])
calc_mean_sd(5000, 6000, cg.mdl$lmp[,1])



# Autocorrelation plots 

# acf computed on the theta counts 

par(mfrow = c(4,2))
acf(fg.mdl$theta[1,1,], lag.max=50, main=expression(paste("Full GS: ", theta)[1][1]), lwd=2, cex = .5)
acf(cg.mdl$theta[1,1,], lag.max=50, main=expression(paste("Collapsed GS: ", theta)[1][1]), lwd=2, cex = .5)

acf(fg.mdl$theta[1,8,], lag.max=50, main=expression(paste("Full GS: ", theta)[8][1]), cex = .5)
acf(cg.mdl$theta[1,8,], lag.max=50, main=expression(paste("Collapsed GS: ", theta)[8][1]), lwd=2, cex = .5)

# acf computed on the beta counts 

acf(fg.mdl$beta[1,1,], lag.max=50, main=expression(paste("Full GS: ", beta)[1][1]), lwd=2, cex = .5)
acf(cg.mdl$beta[1,1,], lag.max=50, main=expression(paste("Collapsed GS: ", beta)[1][1]), lwd=2, cex = .5)

acf(fg.mdl$beta[1,8,], lag.max=50, main=expression(paste("Full GS: ", beta)[1][8]), cex = .5)
acf(cg.mdl$beta[1,8,], lag.max=50, main=expression(paste("Collapsed GS: ", beta)[1][8]), lwd=2, cex = .5)



# Beta plots for the two types of Gibbs samplers 

num.samples <- dim(fg.mdl$beta)[3]
fg.beta <- matrix(0, nrow=K, ncol=V)
cg.beta <- matrix(0, nrow=K, ncol=V)
for (i in (num.samples-100):num.samples){
  fg.beta <- fg.beta + fg.mdl$beta[,,i] / 100;
  cg.beta <- cg.beta + cg.mdl$beta[,,i] / 100;
}

par(mfrow = c(3,1))
color2D.matplot(beta.m, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Synthetic data generating ", beta, " matrix")));
color2D.matplot(fg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Full GS: estimated ", beta, " matrix")));
color2D.matplot(cg.beta, c(0.6, 0), c(0, 0.9), c(0,1), xlab="vocabulary words", ylab="topics", main=expression(paste("Collapsed GS: estimated ", beta, " matrix")));


# Proportion of topic assignments 

total.num.words <- nrow(fg.mdl$Z)
num.samples <- ncol(fg.mdl$Z)
word.id = 2

table(fg.mdl$Z[,word.id])/total.num.words; table(cg.mdl$Z[,word.id])/total.num.words; 



