##############################################################################################
## This script is used to compare the LDA collapsed Gibbs sampler (Griffiths Steyvers 2004)
## and the LDA full Gibbs sampler (Fuentes et al. 2010)
## 
## 
## Last modified on: July 13, 2013 
##############################################################################################

rdata.file     <- "fg_cg_ae33.z.RData"


## Loads packages 

library(ldahp); 

## Initialize variables

set.seed(1983)

K              <- 2 # the number of topics
D              <- 100 # the total number of documents to be generated
V              <- 20 # the vocabulary size
max.iter       <- 40000 # the maximum number of Gibbs iterations
burn.in        <- 2000
spacing        <- 40
lambda.hat     <- 80
gen.eta        <- 3
gen.alpha.v    <- c(3, 3)
gen.eta.v      <- array(gen.eta, c(1, V));                   # symmetric Dirichlet
store.Dir      <- 1                                    # store the \theta and \beta Dirichlet samples ? 

## Generates the synthetic beta.m

beta.m         <- matrix(1e-2, nrow=K, ncol=V)
beta.m[1, ]    <- rdirichlet(1, gen.eta.v);
beta.m[2, ]    <- rdirichlet(1, gen.eta.v);

## Generates documents with a given beta.m

ds             <- generate_docs_fixed_beta(D, lambda.hat, gen.alpha.v, beta.m);



## The full Gibbs sampling

ptm            <- proc.time();
fg.mdl         <- lda_full_c2(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm            <- proc.time();
cg.mdl         <- lda_collapsed_gibbs_c(K, V, ds$wid, ds$doc.N, gen.alpha.v, gen.eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");



##############################################################################################
## The below script is used to compare the proportion of topic assignments (z) 
## from both the collapsed Gibbs sampler (z) and the full Gibbs sampler (z, \beta, \theta) 
##############################################################################################

total.num.words <- nrow(fg.mdl$Z)
num.samples <- ncol(fg.mdl$Z)

fg.cg.ed <- matrix(0, nrow=1, ncol=total.num.words)
for (i in 1:total.num.words){
  ft <- table(fg.mdl$Z[i,])/num.samples; 
  ct <- table(cg.mdl$Z[i,])/num.samples; 
  fg.cg.ed[i] <- sqrt(sum((ft - ct)^2))
}
mean(fg.cg.ed); sd(fg.cg.ed) # mean and s.d. of Euclidean distances between fg and cg topic allocations  


## Saves all objects into a file

save.image(rdata.file)


# word.id = 2
# ft <- table(fg.mdl$Z[word.id,])/num.samples; 
# ct <- table(cg.mdl$Z[word.id,])/num.samples; 
# sqrt(sum((ft - ct)^2))



##############################################################################################
## The below script is used to compare the \theta estimates from 
## the collapsed Gibbs sampler (Griffiths, 2004) and the theta samples 
## from the full Gibbs sampler (Fuentes, 2010) 
##############################################################################################




load("fg_cg_ae108_z.RData")


num.samples <- dim(fg.mdl$theta)[3]
fg.theta <- matrix(0, nrow=K, ncol=D)
cg.theta <- matrix(0, nrow=K, ncol=D)

for (i in 1:num.samples){
  # compute the average of the theta samples which is from the full GS chain 
  fg.theta <- fg.theta + (fg.mdl$theta[,,i] / num.samples);
  
  # compute  the average of the theta estimates which is from the collapsed GS chain  
  cg.theta <- cg.theta + (normalize(cg.mdl$theta[,,i], dim=1) / num.samples);
}

# full GS and collapsed GS theta_d error 
mean(colSums((fg.theta - cg.theta)^2));
sd(colSums((fg.theta - cg.theta)^2));



# two-sample test 
library(ldahp); 
load("fg_cg_ae33_z.RData")

d <- 5
fgt.d <- fg.mdl$theta[, d, ]
cgt.d <- normalize(cg.mdl$theta[, d, ], dim=1) # normalize over the columns 

# t-test using the vector of theta_d vectors 
t.test(fgt.d, cgt.d)

# t-test using the vector of the first element of theta_d vectors 
t.test(fgt.d[1,], cgt.d[1,])

