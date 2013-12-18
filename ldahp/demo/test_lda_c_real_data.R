##############################################################################################
## This script is to test the hyperparameter estimation algorithm  
## using the likelihood ratios on a real dataset. We used the documents 
## extracted from the Wikipedia for this experiment 
##
##############################################################################################


# Sets the working dir
setwd('~/workspace/tm/datasets/whales-tires')
set.seed(1983)
options(digits=2)

library(MCMCpack)
library(ldahp)

## Loads Wikipedia data (from Categories: Whales and Tires) 

vocab <- readLines('whales-tires.ldac.vocab');
documents <- read_docs('whales-tires.ldac');
ds <- vectorize_docs(documents)
doc.N <- calc_doc_lengths(documents)


## Initialize variables
V              <- length(vocab)
K              <- 2                                 # the number of topics
alpha.v        <- c(1, 1)              # symmetric Dirichlet
eta            <- 1 
max.iter       <- 1000                              # the maximum number of Gibbs iterations
burn.in        <- 800
spacing        <- 1
store.Dir      <- 1


## The Gibbs sampling
## Based on the C++ implementation
## Always append the vocabulary ids by 1, because LDA-C consider vocabulary starts at 0

## The full Gibbs sampling

ptm            <- proc.time();
fg.mdl         <- lda_fgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


## The collapsed Gibbs sampling

ptm            <- proc.time();
cg.mdl         <- lda_acgs(K, V, ds$wid+1, doc.N, alpha.v, eta, max.iter, burn.in, spacing, store.Dir);
ptm            <- proc.time() - ptm;
cat("execution time = ", ptm[3], "\n");


#####################################################################################################################
## Testing the LDA library in R 
#####################################################################################################################
library(lda);

setwd('~/workspace/tm/datasets/whales-tires')
set.seed(1983); 

num.topics  <- 2; 
itr         <- 2000; ## Num iterations
av          <- 1;
ev          <- 1;
bp          <- 1600;

documents   <- read.documents(filename = "whales-tires.ldac");
vocab       <- read.vocab(filename = "whales-tires.ldac.vocab")

result      <- lda.collapsed.gibbs.sampler(documents, num.topics, vocab, itr, alpha=av, eta=ev, burnin=bp, trace=2L);
top.topic.words(result$topics, 20, by.score=TRUE);


## Predict new words for the first two documents
predictions <-  predictive.distribution(result$document_sums[,1:2],
                                        result$topics,
                                        0.1, 0.1)

## Use top.topic.words to show the top 5 predictions in each document.
top.topic.words(t(predictions), 5)


attributes(result)

result$assignments



