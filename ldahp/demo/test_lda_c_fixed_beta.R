##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## with a fixed beta matrix (i.e., assume, we have a set of topic 
## distributions before-hand) and we only sample theta in the Gibbs 
## sampling chain. Here, we use the C++ implementation of the Gibbs 
## sampler.  
##################################################################################

library(ldahp);

options(digits=2)
set.seed(1983);
rdata_file          <- "~/Dropbox/lda-hp/results/ldac_fb_cfg02_1.RData" ## to save R data files


## Variables 
## Ref April 24, 2012 mail Dr. Doss
K                   <- 2
D                   <- 100
V                   <- 20
start               <- 0.2
end                 <- 12
interval            <- 0.4
max.iter            <- 9000
burn.in             <- 1000
spacing             <- 100
base.alpha.idx      <- 218  # c(3, 3), the alpha used to run the Gibbs sampler 
gen.alpha.v         <- c(3, 3) # the data generating alpha
lambda.hat          <- 40 # to get the document word count (an auxiliary variable)


alphas              <- gen_meshgrid(start, end, interval) # generates alpha grid (2-D)
base.alpha.v        <- alphas[,base.alpha.idx] # base alpha for the MCMC


## Generates the synthetic beta.m

beta.m              <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, 1:12]     <- array(10,c(1,12))
beta.m[2, 8:20]     <- array(10,c(1,13))


## Generates documents with a given beta.m

ds                  <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m)

## The Gibbs sampling based on the C++ implementation 

ptm                 <- proc.time();
zt.model            <- lda_z_theta_fixed_beta_c(ds, base.alpha.v, beta.m, max.iter, burn.in, spacing, store.Dir=0);
ptm                 <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");

## TODO: there was some memory issues when the objects 
##       are passed back to R from C++. Need to check it out   


thetas <- compute_thetas(ds$did, zt.model$Z, K, D, base.alpha.v)



## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha

log.nu.alphas       <- calc_log_nu_alphas(thetas, alphas)

## Calculates the likelihood ratios

ratios              <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);

## Finds the maximum alpha
max.idx             <- which(ratios == max(ratios)); alphas[,max.idx];
s                   <- sort(ratios, decreasing=T, method="qu", index.return=TRUE); alphas[,s$ix[1:5]];


display_ratios(ratios, start, end, interval);


# ## Saves every object into a file 
# save.image(rdata_file)


# # ###############################################################################################
# # Display results from stored files 
# # ###############################################################################################

postscript(sp_file);
display_ratios(ratios, start, end, interval);
dev.off();

# # ###############################################################################################


###############################################################################################
# Recalculates \theta from the Z values to test 
# if there is any sampling errors  
###############################################################################################
compute.theta.error <- function(did, Z, K, D, base.alpha.v, mc.thetas)
{
  theta.err <- 0; 
  total.N         <- dim(Z)[1];
  sample.count    <- dim(Z)[2];
  thetas          <- array(0, dim=c(K, D, sample.count));
  
  for (iter in 1:sample.count) {
    zid             <- Z[, iter];
    
    # base.alpha.v is used because in the Gibbs sampler we use this 
    theta       <- kronecker(matrix(1, 1, D), base.alpha.v); 
    for (i in 1:total.N) { theta[zid[i], did[i]] <- theta[zid[i], did[i]] + 1; }
    for (d in 1:D) { theta[, d] <- rdirichlet(1, theta[, d]); }
    
    cat("theta.err = ", sum(mc.thetas[,,iter] - theta), "\n");
    theta.err <- theta.err + sum(mc.thetas[,,iter] - theta)
    
    thetas[,,iter] <- theta;
  }
  
  theta.err
}
compute.theta.error(ds$did, zt.model$Z, K, D, base.alpha.v, zt.model$thetas);


log.nu.alphas  <- calc_log_nu_alphas(thetas, alphas)
ratios         <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);
max.idx        <- which(ratios == max(ratios)); alphas[,max.idx];
s              <- sort(ratios, decreasing=T, method="qu", index.return=TRUE); alphas[,s$ix[1:5]];







###############################################################################################


