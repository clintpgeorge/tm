# # --------------------------------------------------------------------------------------
# # Evaluating the Gibbs sampler estimates of theta  
# # --------------------------------------------------------------------------------------


calc.var.original.est <- function(doc.id) {
  ## This function computes the variability 
  ## from document theta, which is used to generate 
  ## the document words, and the estimated 
  ## theta from both the collapsed Gibbs sampler 
  ## and the full Gibbs sampler
  ## 
  ## Arguments: 
  ##  doc.id - document index  
  ## 
  ## Returns: 
  ##  cg.var - the variance between the synthetic theta used to 
  ##  generate the document words and the estimated theta from 
  ##  the collapsed Gibbs sampler (Griffiths 2004)  
  ##
  ##  fg.var - the variance between the synthetic theta used to 
  ##  generate the document words and the estimated theta from 
  ##  the full Gibbs sampler (Fuentes 2010)
  ##
  ##  gen.var - the variance between the synthetic theta used to 
  ##  generate the document words and the theta counts calculated 
  ##  from the synthetically generated document words  
  ## 
  
  # the theta estimate from the collapsed Gibbs sampler 
  cg.theta.est <- cg.mdl$theta[, doc.id] 
  cg.theta.est <- cg.theta.est / sum(cg.theta.est)
  
  # the theta estimate from the Full Gibbs sampler 
  fg.theta.est <- fg.mdl$theta[, doc.id]
  fg.theta.est <- fg.theta.est / sum(fg.theta.est)
  
  # the synthetic theta used to generate the words in the document 
  gen.theta <- ds$theta.samples[, doc.id] 
  
  # the theta counts are counted from the generated words 
  # in the document using the above theta  
  gen.theta.counts <- ds$theta.counts[, doc.id] + gen.alpha.v 
  gen.theta.counts <- gen.theta.counts / sum(gen.theta.counts)   
  
  # Computes variance 
  
  list(cg.var=sqrt(sum((gen.theta - cg.theta.est)^2)), 
       fg.var=sqrt(sum((gen.theta - fg.theta.est)^2)),
       gen.var=sqrt(sum((gen.theta - gen.theta.counts)^2)))
}


# Evaluating estimates 

rdata_file <- "/home/clintpg/results/lda_text_len_d1.RData"
load(rdata_file)

> calc.var.original.est(501)
$cg.var - 0.08736596
$fg.var - 0.08753072
$gen.var - 0.007788693

> calc.var.original.est(1)
$cg.var - 0.1011744
$fg.var - 0.1014368
$gen.var - 0.0759883




