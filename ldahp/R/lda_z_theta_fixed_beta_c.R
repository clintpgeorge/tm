lda_z_theta_fixed_beta_c <- function(ds, alpha.v, beta, max.iter=100, burn.in=0, spacing=1)
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
  
  doc.N        <- ds$doc.N                                            # document word counts
  wid          <- ds$wid                                              # word ids (1 X total.N vector)
  K            <- nrow(beta);                                         # the number of topics
  total.N      <- length(wid);                                        # the total number of word instances 
  n.alpha.v    <- alpha.v / sum(alpha.v);
  zid          <- sample(1:K, total.N, replace=T, prob=n.alpha.v);    # initial random selection of topics for words 
  
  # NOTES: 
  # we substract zid with one because, in C the indexing starts at zero 
  # we assume that the vocab-id also starts at zero
  ret          <- .Call("lda_z_theta_fixed_beta", doc.N, wid-1, zid-1, alpha.v, beta, max.iter, burn.in, spacing, PACKAGE="ldahp");
  
  list(Z=ret$Z+1, thetas=ret$thetas);
  
}