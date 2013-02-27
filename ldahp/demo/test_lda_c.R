##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## and we sample theta and beta in the Gibbs sampling chain. 
## Here, we use the C++ implementation of the Gibbs sampler.  
##################################################################################

library(ldahp);
options(digits=2)
set.seed(1983);

rdata_file     <- "~/Dropbox/lda-hp/results/ldac_cfg01_3.RData"


## Initialize variables
K              <- 2                                    # the number of topics
D              <- 300                                  # the total number of documents to be generated
V              <- 20                                   # the vocabulary size
start          <- 0.2
end            <- 12
interval       <- 0.4
max.iter       <- 11000                                # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 50
lambda.hat     <- 80
gen.alpha.v    <- c(3, 3)                              # symmetric Dirichlet
eta.v          <- array(2, c(1, V));                   # symmetric Dirichlet
base.alpha.idx <- 218                                  # c(3, 3)

alphas         <- gen_meshgrid(start, end, interval)   # generate alpha grid (2-D)
alphas[,base.alpha.idx]                                # base alpha for the MCMC
alpha.v        <- array(alphas[1,base.alpha.idx], dim=c(K, 1));                          # symmetric Dirichlet
eta            <- alphas[2,base.alpha.idx];

set.seed(1983)

## Generates the synthetic beta.m
beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, eta.v);
beta.m[2, ]    <- rdirichlet(1, eta.v);


## Generates documents with a given beta.m
ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);


## The Gibbs sampling

## Based on the C++ implementation
ptm            <- proc.time();
model          <- lda_full_c2(K, V, ds, alpha.v, eta, max.iter, burn.in, spacing);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha
sym.alphas     <- kronecker(matrix(1.0, 1, K), alphas[1,]);
sym.etas       <- kronecker(matrix(1.0, 1, V), alphas[2,]);
log.nu.alphas  <- calc_log_nu_alphas_with_beta(model$thetas, model$betas, sym.alphas, sym.etas);


## Calculates the likelihood ratios
ratios         <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);

## To find the maximum alpha
max.idx        <- which(ratios == max(ratios));
alphas[,max.idx]


s              <- sort(ratios, decreasing = T, method = "qu", index.return=TRUE);
alphas[,s$ix[1:5]];
ratios[s$ix[1:5]];


## Saves every object into a file 
save.image(rdata_file)


## Plots the likelihood ratios (ony for 2-D)
display_ratios(ratios, start, end, interval)


