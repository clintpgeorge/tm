##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## where we sample theta and beta in the Gibbs sampling chain.     
##################################################################################

## includes necessary packages 
library(ldahp);
options(digits=2)
set.seed(1983)


## Initialize variables
rdata_file     <- "lda_fb_cofig01_1.RData"
K              <- 2                                    # the number of topics
D              <- 300                                  # the total number of documents to be generated
V              <- 20                                   # the vocabulary size
start          <- 0.2
end            <- 12
interval       <- 0.4
max.iter       <- 6000                                # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 100
lambda.hat     <- 80
gen.alpha.v    <- c(3, 3)                              # symmetric Dirichlet
eta            <- array(3, c(1, V));                   # symmetric Dirichlet
base.alpha.idx <- 218                                  # c(3, 3)
alphas         <- gen_meshgrid(start, end, interval)   # generate alpha grid (2-D)
base.alpha.v   <- alphas[,base.alpha.idx]              # base alpha for the MCMC

## Generates the synthetic beta.m
beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, eta);
beta.m[2, ]    <- rdirichlet(1, eta);

## Generates documents with a given beta.m
ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);

## The Gibbs sampling
ptm            <- proc.time();

alpha.v        <- array(1, dim=c(K, 1)) * base.alpha.v[1];                          # symmetric Dirichlet
eta            <- base.alpha.v[2];                                                  # symmetric Dirichlet
# model          <- lda_z_theta(K, V, ds, alpha.v, eta, max.iter, burn.in, spacing);
model          <- lda_full(K, V, ds, alpha.v, eta, max.iter, burn.in, spacing);

ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");

## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha
ptm            <- proc.time();

sym.alphas     <- kronecker(matrix(1.0, 1, K), alphas[1,]);
sym.etas       <- kronecker(matrix(1.0, 1, V), alphas[2,]);
log.nu.alphas  <- calc_log_nu_alphas_with_beta(model$thetas, model$betas, sym.alphas, sym.etas);

ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");

## Calculates the likelihood ratios
ratios         <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);


## To find the maximum alpha
max.idx        <- which(ratios == max(ratios));
alphas[,max.idx]
# [1] 4.2 4.2

s              <- sort(ratios, decreasing = T, method = "qu", index.return=TRUE);
alphas[,s$ix[1:5]];
#      [,1] [,2] [,3] [,4] [,5]
# [1,]  4.2  4.6  4.2  4.2  4.6
# [2,]  4.2  3.8  3.8  4.6  4.2

## Saves every object into a file 
save.image(rdata_file)

## Plots the likelihood ratios (ony for 2-D)
# trellis.device(postscript, file="bf-gen=(7,3)-l=80-D=200-G=10000.ps",
#                height=4.8, width=6.5, horiz=F,
#                title="bf-gen=(7,3)-l=80-D=200-G=10000", onefile=T);
display_ratios(ratios, start, end, interval);
# dev.off();



