lda_collapsed_gibbs_fixed_beta <- function(ds, alpha.v, beta, max.iter=100, burn.in=0, spacing=1)
{
  # The LDA collapsed Gibbs sampler  (samples z)
  # This function does not update beta, we just use the given beta to sample z s  
  # 
  # References: 
  #   1. Finding Scienctific Topics Griffiths and Steyvers (2004)
  #   2. David Newman (UCI)'s topic modeling toolkit 
  # 
  # input:
  #   ds       - the corpus       
  #   alpha.v  - hyper parameter vector for theta 
  #   beta     - the beta used to sample the documents   
  #   max.iter - max number of Gibbs iterations to perform 
  #   burn.in  - the burn in period of the Gibbs sampler 
  #   spacing  - spacing between the stored samples (to reduce correlation)
  #
  # return: 
  #   model   - the learned LDA model 
  #
  
  # initializes the variables 
  
  doc.idx <- ds$doc.idx    # document word indices
  doc.N <- ds$doc.N        # document word counts
  wid <- ds$wid            # word ids (1 X total.N vector)
  did <- ds$did            # document ids (1 X total.N vector)
  
  K <- nrow(beta);
  V <- ncol(beta);
  D <- length(doc.N);
  total.N <- length(wid);
  
  valid.samples <- ceiling((max.iter - burn.in) / spacing); # G 
  Z <- matrix(0, nrow=valid.samples, ncol=total.N);
  thetas <- array(0, dim=c(K, D, valid.samples));
  
  theta.counts <- matrix(0, nrow=K, ncol=D);  
  zid <- sample(1:K, total.N, replace=T);                 # initial random selection of topics for the corpus words
  for (i in 1:total.N) { 
    theta.counts[zid[i], did[i]] <- theta.counts[zid[i], did[i]] + 1; 
  }
  
  order <- sample.int(total.N, size=total.N, replace=F);  # for better mixing 
  
  
  # the Gibbs sampler
  count <- 1;
  for (iter in 1:max.iter) {      
    
    for (i in 1:total.N) {
      
      word.idx <- order[i];
      word.id <- wid[word.idx]; 
      doc.id <- did[word.idx];
      z.id <- zid[word.idx];
      mp <- array(0, dim=c(K, 1));
      
      
      theta.counts[z.id, doc.id] <- theta.counts[z.id, doc.id] - 1;
      
      # ref #1: equation 5
      doc.denom <- doc.N[doc.id] - 1 + sum(alpha.v);
      for (k in 1:K){
        mp[k] <- ((theta.counts[k, doc.id] + alpha.v[k]) / doc.denom) * beta[k, word.id]; # here we do not update beta 
      }    
      z.id <- which(rmultinom(1, size=1, prob=mp) == 1); # samples topic index 
      
      theta.counts[z.id, doc.id] <- theta.counts[z.id, doc.id] + 1;
      zid[word.idx] <- z.id;
      
    }    
    
    
    if ((iter > burn.in) && (iter %% spacing == 0)) {
      
      Z[count, ] <- zid;
      
      theta.samples <- matrix(0, nrow=K, ncol=D);
      for (d in 1:D){ theta.samples[, d] <- rdirichlet(1, theta.counts[, d] + alpha.v); }
      thetas[,,count] <- theta.samples;
      
      count <- count + 1;
    }
    
    cat("iter = ", iter, "\n");
    
  }
  
  list(theta.counts=theta.counts, Z=Z, thetas=thetas)
  
}