lda_z_theta_fixed_beta <- function(ds, alpha.v, beta, max.iter=100, burn.in=0, spacing=1)
{
  # The LDA Gibbs sampler: samples z and theta, a given fixed beta    
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
  D <- length(doc.N);                                      # the total number of documents 
  total.N <- length(wid);
  
  zid <- sample(1:K, total.N, replace=T);                  # initial random selection of topics for the corpus words
  
  theta.counts <- matrix(0, nrow=K, ncol=D);
  theta.samples <- matrix(0, nrow=K, ncol=D);
  
  valid.samples <- ceiling((max.iter - burn.in) / spacing); # G 
  Z <- matrix(0, nrow=valid.samples, ncol=total.N);
  thetas <- array(0, dim=c(K, D, valid.samples));
  
  # the Gibbs sampler
  
  count <- 1; 
  for (iter in 1:max.iter){
    
    for (d in 1:D) {
      
      word.idx <- doc.idx[[d]];
      
      ## samples theta 
      Nt <- alpha.v; 
      d.zid <- zid[word.idx];
      for (k in 1:K) { Nt[k] <- Nt[k] + sum(d.zid == k); }
      theta_d <- rdirichlet(1, Nt);
      
      theta.counts[, d] <- Nt; 
      theta.samples[, d] <- theta_d;
      
      
      ## sample the document word instances 
      for (i in 1:doc.N[d]) {
        zid[word.idx[i]] <- which(rmultinom(1, size=1, prob=(beta[, wid[word.idx[i]]] * theta_d)) == 1); 
      }
      
    }  
    
    
    if ((iter > burn.in) && (iter %% spacing == 0)) {
      Z[count, ] <- zid;
      thetas[,,count] <- theta.samples;
      count <- count + 1;
    }
    
    
    cat("iter = ", iter, "\n");
    
  }
  
  list(theta.counts=theta.counts, thetas=thetas, Z=Z);
  
}
