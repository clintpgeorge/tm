##################################################################################
## This is a script to test hyperparameter estimation for the LDA model 
## and we sample theta and beta in the Gibbs sampling chain. 
## Here, we use the C++ implementation of the Gibbs sampler.  
##################################################################################

library(ldahp);
options(digits=2)
set.seed(1983);

rdata_file     <- "~/workspace/tm/data/fg_synth.RData"


## Initialize variables
K              <- 2                                    # the number of topics
D              <- 200                                  # the total number of documents to be generated
V              <- 20                                   # the vocabulary size
start          <- 0.2
end            <- 12
interval       <- 0.4
max.iter       <- 10000                                # the maximum number of Gibbs iterations
burn.in        <- 1000
spacing        <- 10
lambda.hat     <- 80
gen.alpha.v    <- c(3, 3)                              # symmetric Dirichlet
gen.eta.v      <- array(2, c(1, V));                   # symmetric Dirichlet
base.alpha.idx <- 218                                  # c(3, 3)

alphas         <- gen_meshgrid(start, end, interval)   # generate alpha grid (2-D)
alphas[,base.alpha.idx]                                # base alpha for the MCMC
base.alpha.v        <- array(alphas[1,base.alpha.idx], dim=c(K, 1));                          # symmetric Dirichlet
base.eta            <- alphas[2,base.alpha.idx];



## Generates the synthetic beta.m
beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);


## Generates documents with a given beta.m
ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);


## The Gibbs sampling

## Based on the C++ implementation
ptm            <- proc.time();
model          <- lda_full_c2(K, V, ds, base.alpha.v, base.eta, max.iter, burn.in, spacing, 1);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


ret <- compute_thetas_betas(ds$did, ds$wid, model$Z, K, D, V, base.alpha.v, base.eta);


## Calculates nu.alphas (in log scale) using the z and theta samples
## of the MCMC chain ran on the base alpha
sym.alphas     <- kronecker(matrix(1.0, 1, K), alphas[1,]);
sym.etas       <- kronecker(matrix(1.0, 1, V), alphas[2,]);
log.nu.alphas  <- calc_log_nu_alphas_with_beta(ret$thetas, ret$betas, sym.alphas, sym.etas);


## Calculates the likelihood ratios
ratios         <- calc_likelihood_ratios(log.nu.alphas, base.alpha.idx);

## To find the maximum alpha
max.idx        <- which(ratios == max(ratios));
alphas[,max.idx]


s              <- sort(ratios, decreasing = T, method = "qu", index.return=TRUE);
alphas[,s$ix[1:5]];
ratios[s$ix[1:5]];

# 
# ## Saves every object into a file 
# save.image(rdata_file)


## Plots the likelihood ratios (ony for 2-D)
display_ratios(ratios, start, end, interval)




###############################################################################################
# Recalculates \theta and \beta from the Z values and tests 
# if there is any sampling errors  
###############################################################################################

compute.theta.beta.counts.error <- function(did, wid, Z, K, D, V, mc.thetas, mc.betas)
{
  theta.error <- 0; 
  beta.error <- 0; 
  total.N         <- dim(Z)[1];
  sample.count    <- dim(Z)[2];
  
  for (iter in 1:sample.count) {
    zid             <- Z[, iter];
    
    # base.alpha.v and base.eta is used, because the Gibbs sampler uses the same 
    theta       <- matrix(0, nrow=K, ncol=D); 
    beta        <- matrix(0, nrow=K, ncol=V);
    
    for (i in 1:total.N) { 
      beta[zid[i], wid[i]] <- beta[zid[i], wid[i]] + 1; 
    }
    
    for (i in 1:total.N) { 
      theta[zid[i], did[i]] <- theta[zid[i], did[i]] + 1; 
    }
    
    te <- sum(mc.thetas[,,iter] - theta)
    be <- sum(mc.betas[,,iter] - beta)
    cat("theta error = ", te, " beta error = ", be, "\n");

    theta.error <- theta.error + te
    beta.error <- beta.error + be 
    
    
  }
  
}

compute.theta.beta.counts.error(ds$did, ds$wid, model$Z, K, D, V, model$thetas, model$betas)


compute.theta.beta.error <- function(did, wid, Z, K, D, V, base.alpha.v, base.eta, mc.thetas, mc.betas)
{
  theta.error <- 0; 
  beta.error <- 0; 
  total.N         <- dim(Z)[1];
  sample.count    <- dim(Z)[2];
  
  for (iter in 1:sample.count) {
    zid             <- Z[, iter];
    
    # base.alpha.v and base.eta is used, because the Gibbs sampler uses the same 
    theta       <- kronecker(matrix(1, 1, D), base.alpha.v); 
    beta        <- matrix(base.eta, nrow=K, ncol=V);
    
    for (i in 1:total.N) { 
      beta[zid[i], wid[i]] <- beta[zid[i], wid[i]] + 1; 
    }
    for (k in 1:K) { 
      beta[k,] <- rdirichlet(1, beta[k,]); 
    } 
    
    for (i in 1:total.N) { 
      theta[zid[i], did[i]] <- theta[zid[i], did[i]] + 1; 
    }
    for (d in 1:D) { 
      theta[, d] <- rdirichlet(1, theta[, d]); 
    }
    
    te <- sum(mc.thetas[,,iter] - theta)
    be <- sum(mc.betas[,,iter] - beta)
    cat("theta error = ", te, " beta error = ", be, "\n");
    
    theta.error <- theta.error + te
    beta.error <- beta.error + be 
    
    
  }
  
  cat("avg theta error = ", theta.error/sample.count, " avg beta error = ", beta.error/sample.count, "\n");
}

compute.theta.beta.error(ds$did, ds$wid, model$Z, K, D, V, base.alpha.v, base.eta, model$thetas, model$betas)
