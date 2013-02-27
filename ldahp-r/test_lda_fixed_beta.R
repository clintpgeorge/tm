##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## with a fixed beta matrix (i.e., assume, we have a set of topic 
## distributions before-hand) and we only sample theta in the Gibbs 
## sampling chain.     
##################################################################################

## includes necessary packages 
library(ldahp)
options(digits=2)
set.seed(1983)


## Initializes variables

rdata_file <- "lda_fb_cofig01_1.RData"

##  Ref April 24, 2012 mail Dr. Doss
K <- 2 # the number of topics
D <- 10 # the total number of documents to be generated
V <- 20 # the vocabulary size
start <- 0.2
end <- 12
interval <- 0.4
max.iter <- 6000 # the maximum number of Gibbs iterations
burn.in <- 1000
spacing <- 100
gen.alpha.v <- c(3, 3)
lambda.hat <- 80
base.alpha.idx <- 218 # (3, 3)

## generates alpha grid (2-D)
alphas <- gen_meshgrid(start, end, interval) 
base.alpha.v <- alphas[,base.alpha.idx] # base alpha for the MCMC


## Generates the synthetic beta.m
beta.m <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, 1:12] <- array(10, c(1,12))
beta.m[2, 8:20] <- array(10, c(1,13))


## Generates documents with a given beta.m
ds <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);

## The Gibbs sampling
ptm <- proc.time();
zt.model <- lda_z_theta_fixed_beta(ds, base.alpha.v, beta.m, max.iter, burn.in, spacing);
# cb.model <- lda_collapsed_gibbs_fixed_beta(ds, base.alpha.v, beta.m, max.iter, burn.in, spacing);
ptm <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha
log.nu.alphas <- calc_log_nu_alphas(zt.model$thetas, alphas);

## Calculates the likelihood ratios
ratios <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);


## Finds the maximum alpha
max.idx <- which(ratios == max(ratios));
alphas[,max.idx]

s <- sort(ratios, decreasing=T, method="qu", index.return=TRUE);
alphas[,s$ix[1:5]];

## Saves every object into a file 
save.image(rdata_file)



## Plots the likelihood ratios (ony for 2-D)
display_ratios(ratios, start, end, interval);



