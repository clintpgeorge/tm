lda_collapsed_gibbs_c <- function(K, V, ds, alpha.v, eta, max.iter=100, burn.in=0, spacing=1)
{
  # The LDA Gibbs sampler  (samples z)
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
  
  
  doc.N        <- ds$doc.N;                                         # document word counts
  wid          <- ds$wid;                                           # word ids (1 X total.N vector)
  total.N      <- length(wid);                                      # the total number of word instances 
  n.alpha.v    <- alpha.v / sum(alpha.v);
  zid          <- sample(1:K, total.N, replace=T, prob=n.alpha.v);  # initial selection of topics for words
  
  # NOTES: 
  # we substract zid with one because, in C the indexing starts at zero 
  # we assume that the vocab-id also starts at zero
  ret          <- .Call("lda_collapsed_gibbs", K, V, doc.N, wid-1, zid-1, alpha.v, eta, max.iter, burn.in, spacing, PACKAGE="ldahp");
  
  list(Z=ret$Z+1, thetas=ret$thetas, betas=ret$betas, lmp=ret$lmp);
  
}