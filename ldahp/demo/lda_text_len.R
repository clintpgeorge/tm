library(plotrix);
library(ldahp);


generate_lda_docs <- function(doc.N, alpha.v, beta)
{
  
  # Generates document words using the LDA generative process given a beta. 
  # It's used to test the correctness of the Gibbs sampling algrithms.
  # 
  # Inputs: 
  #   doc.N      - the document counts  
  #   alpha.v    - the vector of Dirichlet hyperparameters (K X 1) for document topic mixtures 
  #   beta       - the beta matrix (counts) for topic word probabilities (K x V format)
  # 
  
  K <- nrow(beta);                                 # the number of topics 
  V <- ncol(beta)                                  # the vocabulary size 
  D <- length(doc.N);
  theta.counts <- matrix(0, nrow=K, ncol=D);       # stores document topic word counts  
  beta.counts <- matrix(0, nrow=K, ncol=V);        # stores topic word counts 
  theta.samples <- matrix(0, nrow=K, ncol=D);      
  
  did <- c();
  wid <- c();
  zid <- c();
  
  
  num <- 1;
  doc.idx <- vector("list", D);
  
  for (d in 1:D)
  {
    ptm <- proc.time();
    
    theta.samples[,d] <- rdirichlet(1, alpha.v);
    
    did <- cbind(did, array(1, c(1, doc.N[d])) * d); # document instances 
    
    z_d <- c(); 
    indices <- c();
    for (i in 1:doc.N[d]){
      
      z_dn <- which(rmultinom(1, size=1, prob=theta.samples[,d]) == 1); # samples topic 
      w_dn <- which(rmultinom(1, size=1, beta[z_dn,]) == 1); # samples word
      
      wid <- cbind(wid, w_dn); 
      z_d <- cbind(z_d, z_dn);  
      indices <- cbind(indices, num);
      
      num <- num + 1;            
      
    }
    doc.idx[[d]] <- as.integer(indices); # stores the document word indices     
    
    
    theta.counts[, d] <- calc_topic_counts(z_d, K); # calculates the document topic counts
    
    zid <- cbind(zid, z_d); 
    
    ptm <- proc.time() - ptm;
    cat("document = ", d, " time = ", ptm[3], " # words = ", doc.N[d], "\n");
  }
  
  
  total.N <- sum(doc.N);
  
  for (i in 1:total.N){ 
    beta.counts[zid[i], wid[i]] <- beta.counts[zid[i], wid[i]] + 1;
  }
  
  list(did=as.vector(did), wid=as.vector(wid), zid=as.vector(zid), 
       theta.counts=theta.counts, beta.counts=beta.counts, theta.samples=theta.samples, 
       total.N=total.N, doc.N=doc.N, doc.idx=doc.idx);
  
}


# # --------------------------------------------------------------------------------------
# # Generates synthetic data and Gibbs sampler runs   
# # --------------------------------------------------------------------------------------


## Initialize variables

set.seed(1983)

K              <- 30 # the number of topics
D              <- 501 # the total number of documents to be generated
V              <- 50 # the vocabulary size
max.iter       <- 20000 # the maximum number of Gibbs iterations
burn.in        <- 5000
spacing        <- 10
lambda.hat     <- 50
gen.eta        <- 3
gen.alpha      <- 3
gen.alpha.v    <- array(gen.alpha, c(1, K)); 
gen.eta.v      <- array(gen.eta, c(1, V));                   # symmetric Dirichlet
rdata_file     <- "lda_text_len_d2.RData"


## Generates the synthetic beta.m
beta.m         <- matrix(1e-2, nrow=K, ncol=V)
for (i in 1:K){ beta.m[i, ]    <- rdirichlet(1, gen.eta.v); }


## Generates documents with a given beta.m
doc.N          <- array(lambda.hat, dim=c(D, 1));
doc.N[D]       <- 10000
ds             <- generate_lda_docs(doc.N, gen.alpha.v, beta.m);



## The full Gibbs sampling

ptm            <- proc.time();
fg.mdl         <- lda_fg(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


# ## The collapsed Gibbs sampling

ptm               <- proc.time();
cg.mdl            <- lda_cg(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing);
ptm               <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


save.image(rdata_file)





