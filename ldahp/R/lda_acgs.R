lda_acgs <- function(K, V, wid, doc.N, alpha.v, eta, max.iter=100, burn.in=0, spacing=1, store.Dir=1)
{
  # The LDA Augmented Collapsed Gibbs sampler  (samples Z, Beta, and Theta) and 
  # the Collapsed Gibbs sampler (sampler Z)
  # 
  # input:
  #   K        - the number of topics in the corpus (Guessed) 
  #   V        - the vocabulary size 
  #   wid      - the vocabulary ids of every word instance in each corpus document  
  #              (1 X total.N vector). We assume vocabulary id starts with 1      
  #   doc.N    - the document lengths   
  #   alpha.v  - the hyper parameter vector for theta 
  #   eta      - beta matrix smoothing parameter 
  #   max.iter - the max number of Gibbs iterations to be performed  
  #   burn.in  - the burn in period of the Gibbs sampler 
  #   spacing  - the spacing between the stored samples (to reduce correlation)
  #   store.Dir- if 1 ACGS is performed else CGS is performed    
  #
  # return: 
  #   model    - the learned LDA model 
  #
  
  # initializes the variables 
  
  total.N      <- length(wid);                                      # the total number of word instances 
  n.alpha.v    <- alpha.v / sum(alpha.v);
  zid          <- sample(1:K, total.N, replace=T, prob=n.alpha.v);  # initial selection of topics for words
  
  # NOTES: 
  # we substract zid with one because, in C the indexing starts at zero 
  # we assume that the vocab-id also starts at zero
  ret          <- .Call("lda_acgs", K, V, wid-1, doc.N, zid-1, alpha.v, eta, max.iter, burn.in, spacing, store.Dir, PACKAGE="ldahp");
  
  list(Z=ret$Z+1, theta=ret$thetas, beta=ret$betas, lmp=ret$lmp, lp=ret$lp);
  
}