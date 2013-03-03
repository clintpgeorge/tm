lda_full_c <- function(K, V, ds, alpha.v, eta, max.iter=100, burn.in=0, spacing=1, store.Dir=0)
{
    # The LDA Gibbs sampler   (samples z, beta, and theta)
    # 
    # input:
    #   K        - the number of topics 
    #   V        - the vocabulary size 
    #   ds       - the corpus (read using read.ldac.documents from lda.R)
    #   alpha.v  - hyper parameter vector for theta 
    #   eta      - beta matrix smoothing parameter 
    #   max.iter - max number of Gibbs iterations to perform 
    #   burn.in  - the burn in period of the Gibbs sampler 
    #   spacing  - spacing between the stored samples (to reduce correlation)
    #   store.Dir- if 0 the sampler does not save theta and beta samples  
    #
    # return: 
    #   model     - the learned LDA model 
    #
    
    # initializes the variables 
    
    total.N      <- sum(ds$doc.lengths);                                      # the total number of word instances 
    n.alpha.v    <- alpha.v / sum(alpha.v);
    zid          <- sample(1:K, total.N, replace=T, prob=n.alpha.v);          # initial selection of topics for words

    # NOTES: 
    # we substract zid with one because, in C the indexing starts at 0 
    # we assume that the vocab-id also starts at zero. Here the corpus
    # is assumed to be read from the LDAC formatted file using 
    # read.ldac.documents function from lda.R
    ret          <- .Call("lda_full", K, V, ds$doc.N, ds$docs, zid-1, alpha.v, eta, max.iter, burn.in, spacing, store.Dir, PACKAGE="ldahp");
    
    list(Z=ret$Z+1, thetas=ret$thetas, betas=ret$betas, lmp=ret$lmp);

}




