################################################################################################
## optimizing the hyperparameters using optimx package
################################################################################################

library(optimx);
library(ldahp);
options(digits=2)


## CASE I: FIXED BETA 

alpha <- c(1, 1); # initial alpha
mc.thetas <- zt.model$thetas; # from the lda gibbs sampler  
base.alpha <- c(3, 3); # the base alpha used by the Gibbs sampler 

# To test the function 
# exp(-neg_log_partial_likelihood_ratio(alpha, mc.thetas, base.alpha));

ans1 <- optimx(par=alpha, fn=neg_log_partial_likelihood_ratio, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
               method="CG", itnmax=100, hessian=FALSE,control=list(trace=1), mc.thetas, base.alpha);
print(ans1);



## CASE I: SAMPLES BETA 

alpha_eta <- c(1, 1, 1); # initial alpha
mc.thetas <- model$thetas; # from the lda gibbs sampler  
mc.betas <- model$betas; # from the lda gibbs sampler
base.alpha_eta <- c(1, 1, 1); # the base alpha used by the Gibbs sampler 

exp(-neg_log_likelihood_ratio(alpha_eta, mc.thetas, mc.betas, base.alpha_eta));

ans1 <- optimx(par=alpha_eta, fn=neg_log_likelihood_ratio, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
               method="CG", itnmax=100, hessian=FALSE, control=list(trace=1), mc.thetas, mc.betas, base.alpha_eta);
print(ans1);





## CASE I: SAMPLES BETA (SYMMETRIC HYPER)


alpha_eta <- c(1,1); # initial alpha
mc.thetas <- model$thetas; # from the lda gibbs sampler  
mc.betas <- model$betas; # from the lda gibbs sampler
base.alpha_eta <- c(3,3); # the base alpha used by the Gibbs sampler 

exp(-neg_log_likelihood_ratio_symmetric(alpha_eta, mc.thetas, mc.betas, base.alpha_eta));

ans1 <- optimx(par=alpha_eta, fn=neg_log_likelihood_ratio_symmetric, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
               method="CG", itnmax=100, hessian=FALSE,control=list(trace=1), mc.thetas, mc.betas, base.alpha_eta);
print(ans1);
