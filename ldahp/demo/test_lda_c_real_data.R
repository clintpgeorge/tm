##############################################################################################
## This script is to test the hyperparameter estimation algorithm  
## using the likelihood ratios on a real dataset. We used the documents 
## extracted from the Wikipedia for this experiment 
##
##############################################################################################


# Sets the working dir
setwd("~/Dropbox/lda-hp/r") 

library(MCMCpack)
library(lda.hp);
library(optimx);

source('lda.R')
options(digits=2)


## Initialize variables
# wp 11000: 1830s, \alpha = c(0.36, 0.54, 1.00), b = 1.00 

vocab          <- readLines('~/Dropbox/lda-data/wp/wp.ldac.vocab');
# ds             <- read.ldac.documents2('~/Dropbox/lda-data/20120604/ap.ldac');
# vocab          <- readLines('~/Dropbox/lda-data/20120604/ap.ldac.vocab');


V              <- length(vocab)
K              <- 10                                 # the number of topics
alpha.v        <- array(1, dim=c(K, 1))              # symmetric Dirichlet
eta            <- 1 
max.iter       <- 100                                # the maximum number of Gibbs iterations
burn.in        <- 50
spacing        <- 1

## The Gibbs sampling

set.seed(1983)

## Based on the C++ implementation
ds             <- read.ldac.documents('~/Dropbox/lda-data/wp/wp.ldac');
model          <- lda_full_c(K, V, ds, alpha.v, eta, max.iter, burn.in, spacing);


ptm            <- proc.time();

ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


ds2             <- read.ldac.documents.Gibbs('~/Dropbox/lda-data/wp/wp.ldac');
ptm            <- proc.time();
model2          <- lda_full_c2(K, V, ds2, alpha.v, eta, max.iter, burn.in, spacing);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


#####################################################################################################################
library(lda);
setwd('~/Dropbox/lda-data/20120604/')

set.seed(1983); 

num.topics  <- 10; 
itr         <- 2000; ## Num iterations
av          <- 1;
ev          <- 1;
bp          <- 1600;

documents   <- read.documents(filename = "ap.ldac");
vocab       <- read.vocab(filename = "ap.ldac.vocab")

result      <- lda.collapsed.gibbs.sampler(documents, num.topics, vocab, itr, alpha=av, eta=ev, burnin=bp, trace=2L);
top.topic.words(result$topics, 10, by.score=TRUE);


## Predict new words for the first two documents
predictions <-  predictive.distribution(result$document_sums[,1:2],
                                        result$topics,
                                        0.1, 0.1)

## Use top.topic.words to show the top 5 predictions in each document.
top.topic.words(t(predictions), 5)


attributes(ds)

