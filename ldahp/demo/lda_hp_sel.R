##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## and we sample theta and beta in the Gibbs sampling chain. 
## Here, we use the C++ implementation of the Gibbs sampler.  
## This can support heyperparametrs taken from a 2-dimensional plane.
##################################################################################
rm(list=ls());

# load the necessary libraries
library(ldahp);
options(digits=2)
set.seed(1983);


## Initialize variables


## h = (eta, alpha) = (3, 3)
rdata.file     <- "fg_ea33.RData"
gen.eta        <- 3 # symmetric Dirichlet
gen.alpha      <- 3 # symmetric Dirichlet
base.alpha.idx <- 855 


# ## h = (eta, alpha) = (7, 7)
# rdata.file     <- "fg_ea77.RData"
# gen.eta        <- 7 # symmetric Dirichlet
# gen.alpha      <- 7 # symmetric Dirichlet
# base.alpha.idx <- 2075 
# 
# ## h = (eta, alpha) = (7, 3)
# rdata.file     <- "fg_ea73.RData"
# gen.eta        <- 7 # symmetric Dirichlet
# gen.alpha      <- 3 # symmetric Dirichlet
# base.alpha.idx <- 875 
# 
# ## h = (eta, alpha) = (10, 8)
# rdata.file     <- "fg_ea108.RData"
# gen.eta        <- 10 # symmetric Dirichlet
# gen.alpha      <- 8 # symmetric Dirichlet
# base.alpha.idx <- 2390 


K              <- 2 # the number of topics
D              <- 1000 # the total number of documents to be generated
V              <- 20 # the vocabulary size
start          <- 0.2
end            <- 12
interval       <- 0.2
max.iter       <- 61000 # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 50
lambda.hat     <- 80
store.Dir      <- 1

gen.alpha.v    <- array(gen.alpha, c(K, 1));          
gen.eta.v      <- array(gen.eta, c(1, V));                   

alphas         <- gen_meshgrid(start, end, interval)   # generate alpha grid (2-D)
alphas[,base.alpha.idx]                                # base alpha for the MCMC

base.alpha.v   <- array(alphas[1,base.alpha.idx], dim=c(K, 1));                          # symmetric Dirichlet
base.eta       <- alphas[2,base.alpha.idx];



## Generates the synthetic beta.m
beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);


## Generates documents with a given beta.m
ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);


## The Gibbs sampling

## Based on the C++ implementation
ptm            <- proc.time();
model          <- lda_full_c2(K, V, ds$wid, ds$doc.N, base.alpha.v, base.eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


# ptm            <- proc.time();
# ret <- compute_thetas_betas(ds$did, ds$wid, model$Z, K, D, V, base.alpha.v, base.eta);
# ptm            <- proc.time() - ptm;
# cat("execution time = ", ptm[3], "\n");

## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha
sym.alphas     <- kronecker(matrix(1.0, 1, K), alphas[1,]);
sym.etas       <- kronecker(matrix(1.0, 1, V), alphas[2,]);
# log.nu.alphas  <- calc_log_nu_alphas_with_beta(ret$thetas, ret$betas, sym.alphas, sym.etas);
log.nu.alphas  <- calc_log_nu_alphas_with_beta(model$thetas, model$betas, sym.alphas, sym.etas);

## Calculates the likelihood ratios
ratios         <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);

## To find the maximum alpha

s              <- sort(ratios, decreasing = T, method = "qu", index.return=TRUE);
alphas[,s$ix[1:20]];
ratios[s$ix[1:20]];

## Saves every object into a file 

save.image(rdata.file)

