################################################################################################
## The following functions are used for the hyperparameter optimization using optimx package
## Added on: May 31, 2012 
################################################################################################



neg_log_partial_likelihood_ratio <- function(alpha, mc.thetas, base.alpha){
  # Computes the LDA negative log likelihood ratio from a fixed beta  
  # (no topic sampling) Gibbs sampler output. This function is used as 
  # an argument for the R optimx package.  
  # 
  # Arguments:
  #   alpha - the hyperparameter vector (1 x K vector) to be optimized  
  #   mc.thetas - this contains sampled thetas (K x D x G matrix) from the LDA Markov chain   
  #   base.alpha - the hyperparameter vector (1 x K vector) used by the LDA Gibbs sampler 
  # 
  # Returns: 
  #   likelihood ratio (float) 
  #    
  
  G <- dim(mc.thetas)[3]; # number of samples from MC      
  ln.const <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - lgamma(sum(base.alpha)) + sum(lgamma(base.alpha)); # these are const for all i's 
  alpha.diff <- alpha - base.alpha
  
  
  dr <- 0;
  for (i in 1:G){
    dr <- dr + exp( sum( ln.const + (alpha.diff  %*% log(mc.thetas[,,i])) ) );
  }
  
  -log(dr / G); # to make it as a minimization problem 
  
}


neg_log_likelihood_ratio <- function(alpha_eta, mc.thetas, mc.betas, base.alpha_eta){
  # Computes the LDA negative log likelihood ratio from an LDA full
  # (no topic sampling) Gibbs sampler output. This function is used as 
  # an argument for the R optimx package.  
  # 
  # Arguments:
  #   alpha_eta - the hyperparameter vector (1 x K+1 vector) to be optimized (the last element is eta)  
  #   mc.thetas - this contains sampled thetas (K x D x G matrix) from the LDA Markov chain   
  #   mc.betas - this contains sampled betas (K x V x G matrix) from the LDA Markov chain   
  #   base.alpha_eta - the hyperparameter vector (1 x K+1 vector) used by the LDA Gibbs sampler (the last element is eta)  
  # 
  # Returns: 
  #   likelihood ratio (float) 
  #
  
  len_ae <- length(alpha_eta);
  alpha <- alpha_eta[-len_ae];
  eta <- alpha_eta[len_ae];
  base.alpha <- base.alpha_eta[-len_ae];
  base.eta <- base.alpha_eta[len_ae];
  G <- dim(mc.thetas)[3]; # number of samples from MC 
  V <- dim(mc.betas)[2]; # number of vocabulary terms 
  eta.v <- array(eta, dim=c(V, 1));
  base.eta.v <- array(base.eta, dim=c(V, 1));
  
  lnc.alpha <- lgamma(sum(alpha)) - sum(lgamma(alpha)) - lgamma(sum(base.alpha)) + sum(lgamma(base.alpha));
  lnc.eta <- lgamma(V * eta) - V * lgamma(eta) - lgamma(V * base.eta) + V * lgamma(base.eta);
  
  dr <- 0;
  for (i in 1:G){
    dr <- dr + exp( sum( ((alpha - base.alpha) %*% log(mc.thetas[,,i])) + lnc.alpha ) + sum( (log(mc.betas[,,i]) %*% (eta.v - base.eta.v)) + lnc.eta ) );
  }
  
  -log(dr / G); # to make it as a minimization problem 
  
}

neg_log_likelihood_ratio_symmetric <- function(alpha_eta, mc.theta, mc.beta, base.alpha_eta){
  # Computes the LDA negative log likelihood ratio from an LDA full
  # (no topic sampling) Gibbs sampler output, with symmetric hyperparameters. 
  # This function is used as an argument for the R optimx package.  
  # 
  # Arguments:
  #   alpha_eta - the symmetric hyperparameter vector (1 x 2 vector) to be optimized (the last element is eta)  
  #   mc.thetas - this contains sampled thetas (K x D x G matrix) from the LDA Markov chain   
  #   mc.betas - this contains sampled betas (K x V x G matrix) from the LDA Markov chain   
  #   base.alpha_eta - the symmetric hyperparameter vector (1 x 2 vector) used by the LDA Gibbs sampler (the last element is eta)  
  # 
  # Returns: 
  #   likelihood ratio (float) 
  #
  
  alpha <- alpha_eta[1];
  eta <- alpha_eta[2];
  base.alpha <- base.alpha_eta[1];
  base.eta <- base.alpha_eta[2];
  G <- dim(mc.thetas)[3]; # number of samples from MC 
  K <- dim(mc.thetas)[1]; # number of topics in the corpus 
  V <- dim(mc.betas)[2]; # number of vocabulary terms 
  
  alpha.v <- array(alpha, dim=c(1, K));
  base.alpha.v <- array(base.alpha, dim=c(1, K));
  eta.v <- array(eta, dim=c(V, 1));
  base.eta.v <- array(base.eta, dim=c(V, 1));
  
  lnc.alpha <- lgamma(K * alpha) - K * lgamma(alpha) - lgamma(K * base.alpha) + K * lgamma(base.alpha);
  lnc.eta <- lgamma(V * eta) - V * lgamma(eta) - lgamma(V * base.eta) + V * lgamma(base.eta);
  
  dr <- 0;
  for (i in 1:G){
    dr <- dr + exp( sum( ((alpha - base.alpha) %*% log(mc.thetas[,,i])) + lnc.alpha ) + sum( (log(mc.betas[,,i]) %*% (eta.v - base.eta.v)) + lnc.eta ) );
  }
  
  -log(dr / iter);
  
}
