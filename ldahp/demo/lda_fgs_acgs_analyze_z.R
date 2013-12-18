##############################################################################################
## This script is used to compare the Augmented Collapsed Gibbs sampler
## and the full Gibbs sampler of LDA based on the z sample (word topic 
## allocations) 
## 
## 
## Last modified on: Sept 27, 2013 
##############################################################################################

rdata.file <- "fg_cg_ea73_pvalues.RData"

## Loads packages 

library(ldahp); 

## Initialize variables

set.seed(1983)

K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 50000 # the maximum number of Gibbs iterations
burn.in        <- 10000
spacing        <- 40
lambda.hat     <- 80
gen.eta        <- 7
gen.alpha      <- 3
gen.alpha.v    <- array(gen.alpha, c(1, K));           # symmetric Dirichlet
gen.eta.v      <- array(gen.eta, c(1, V));             # symmetric Dirichlet
store.Dir      <- 1                                    # store the \theta and \beta Dirichlet samples ? 

## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);



## The full Gibbs sampling

ptm            <- proc.time();
fg.mdl         <- lda_fgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm            <- proc.time();
cg.mdl         <- lda_acgs(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");

##############################################################################################
## 2-samples t-test on z1 
##############################################################################################

total.num.words <- nrow(fg.mdl$Z)
num.samples <- ncol(fg.mdl$Z)

p.values <- matrix(0, nrow=1, ncol=total.num.words)
for (i in 1:total.num.words){
  fgs.z1 <- fg.mdl$Z[i,]
  acgs.z1 <- cg.mdl$Z[i,]
  
  fgs.z1[fgs.z1 == 2] <- 0
  acgs.z1[acgs.z1 == 2] <- 0
  
  t.res <- t.test(fgs.z1, acgs.z1)
  p.values[i] <- t.res$p.value
}

## Saves all objects into a file

save.image(rdata.file)





# ##############################################################################################
# ## The below script is used to compare the proportion of topic assignments (z) 
# ## from both the collapsed Gibbs sampler (z) and the full Gibbs sampler (z, \beta, \theta) 
# ##############################################################################################
# 
# total.num.words <- nrow(fg.mdl$Z)
# num.samples <- ncol(fg.mdl$Z)
# 
# fg.cg.ed <- matrix(0, nrow=1, ncol=total.num.words)
# for (i in 1:total.num.words){
#   ft <- table(fg.mdl$Z[i,])/num.samples; 
#   ct <- table(cg.mdl$Z[i,])/num.samples; 
#   fg.cg.ed[i] <- sqrt(sum((ft - ct)^2))
# }
# mean(fg.cg.ed); sd(fg.cg.ed) # mean and s.d. of Euclidean distances between fg and cg topic allocations  
# 
# 
# 
# 
# # word.id = 2
# # ft <- table(fg.mdl$Z[word.id,])/num.samples; 
# # ct <- table(cg.mdl$Z[word.id,])/num.samples; 
# # sqrt(sum((ft - ct)^2))
# 
# 
# 
# ##############################################################################################
# ## The below script is used to compare the \theta estimates from 
# ## the collapsed Gibbs sampler (Griffiths, 2004) and the theta samples 
# ## from the full Gibbs sampler (Fuentes, 2010) 
# ##############################################################################################
# 
# 
# 
# 
# load("fg_cg_ae108_z.RData")
# 
# 
# num.samples <- dim(fg.mdl$theta)[3]
# fg.theta <- matrix(0, nrow=K, ncol=D)
# cg.theta <- matrix(0, nrow=K, ncol=D)
# 
# for (i in 1:num.samples){
#   # compute the average of the theta samples which is from the full GS chain 
#   fg.theta <- fg.theta + (fg.mdl$theta[,,i] / num.samples);
#   
#   # compute  the average of the theta estimates which is from the collapsed GS chain  
#   cg.theta <- cg.theta + (normalize(cg.mdl$theta[,,i], dim=1) / num.samples);
# }
# 
# # full GS and collapsed GS theta_d error 
# mean(sqrt(colSums((fg.theta - cg.theta)^2)));
# sd(sqrt(colSums((fg.theta - cg.theta)^2)));
# 
# 
# 
