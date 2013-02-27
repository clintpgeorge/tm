lda_full <- function(K, V, ds, alpha.v, eta, max.iter=100, burn.in=0, spacing=1)
{
  # The LDA Gibbs sampler  (samples z and theta)
  # 
  # input:
  #   K        - the number of topics 
  #   V        - the vocabulary size 
  #   ds       - the corpus     
  #   alpha.v  - hyper parameter vector for theta 
  #   eta      - beta matrix smoothing parameter 
  #   max.iter - max number of Gibbs iterations to perform 
  #   burn.in  - the burn in period of the Gibbs sampler 
  #   spacing  - spacing between the stored samples (to reduce correlation)
  #
  # return: 
  #   model     - the learned LDA model 
  #
  
  # initializes the variables 
  
  
  doc.idx <- ds$doc.idx    # document word indices
  doc.N <- ds$doc.N        # document word counts
  wid <- ds$wid            # word ids (1 X total.N vector)
  did <- ds$did            # document ids (1 X total.N vector)
  
  D <- length(doc.N);      # the total number of documents
  total.N <- length(wid);  # the total number of word instances 
  
  valid.samples <- ceiling((max.iter - burn.in) / spacing); # G 
  
  Z <- matrix(0, nrow=valid.samples, ncol=total.N);
  thetas <- array(0, dim=c(K, D, valid.samples));
  betas <- array(0, dim=c(K, V, valid.samples));
  
  n.alpha.v <- alpha.v / sum(alpha.v);
  zid <- sample(1:K, total.N, replace=T, prob=n.alpha.v);                 # initial random selection of topics for the corpus words
  
  beta.counts <- matrix(0, nrow=K, ncol=V);
  beta <- matrix(0, nrow=K, ncol=V);
  for (i in 1:total.N) { 
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1; 
  }
  
  # the Gibbs sampler
  count <- 1
  for (iter in 1:max.iter) {      
    
    theta.samples <- matrix(0, nrow=K, ncol=D);
    for (d in 1:D) {
      
      word.idx <- doc.idx[[d]];
      
      # samples beta 
      for (k in 1:K) { beta[k,] <- rdirichlet(1, beta.counts[k,] + eta); }            
      
      
      ## samples theta 
      Nt <- array(0, dim=c(K, 1)); 
      d.zid <- zid[word.idx];
      for (k in 1:K) { Nt[k] <- Nt[k] + sum(d.zid == k); }
      theta_d <- rdirichlet(1, Nt + alpha.v);
      theta.samples[, d] <- theta_d;
      
      
      # excludes document d's word-topic counts 
      for (i in 1:doc.N[d]){
        beta.counts[zid[word.idx[i]], wid[word.idx[i]]] <- beta.counts[zid[word.idx[i]], wid[word.idx[i]]] - 1; 
      }
      
      ## samples word instances 
      for (i in 1:doc.N[d]) {
        zid[word.idx[i]] <- which(rmultinom(1, size=1, prob=(beta[, wid[word.idx[i]]] * theta_d)) == 1); 
      }
      
      # updates document d's word-topic counts 
      for (i in 1:doc.N[d]){
        beta.counts[zid[word.idx[i]], wid[word.idx[i]]] <- beta.counts[zid[word.idx[i]], wid[word.idx[i]]] + 1; 
      }
      
    }    
    
    
    if ((iter > burn.in) && (iter %% spacing == 0)) {
      Z[count, ] <- zid;
      thetas[,,count] <- theta.samples;
      betas[,,count] <- beta;
      count <- count + 1;
    }
    
    
    cat("iter = ", iter, "\n");
    
  }
  
  list(Z=Z, thetas=thetas, betas=betas)
  
}