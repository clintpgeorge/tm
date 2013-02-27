calc_log_nu_alphas <- function (thetas, alphas) {
  
  sample.count   <- dim(thetas)[3];
  D          <- dim(thetas)[2];
  num.grids  <- ncol(alphas);
  nu.alphas  <- matrix(0, nrow=num.grids, ncol=sample.count);
  
  for (ng in 1:num.grids){
    
    
    alpha.v  <- alphas[, ng];        
    ln.nu.alpha.C <- D * (lgamma(sum(alpha.v)) - sum(lgamma(alpha.v))); # calculates the constant part of the likilihood ratio
    
    for (iter in 1:sample.count){
      nu.alphas[ng, iter] <- sum(t(alpha.v) %*% log(thetas[,,iter])) + ln.nu.alpha.C;  
    }
    
    cat("grid = ", ng, "\n");
    
  } 
  
  nu.alphas;    
  
}

calc_log_nu_alphas_with_beta <- function (thetas, betas, alphas, etas) {
  
  sample.count   <- dim(thetas)[3];
  D          <- dim(thetas)[2];
  num.grids  <- nrow(alphas); # G x K matrix 
  K          <- ncol(alphas);
  V          <- ncol(etas); # G x V matrix 
  nu.alphas  <- matrix(0, nrow=num.grids, ncol=sample.count);
  
  for (ng in 1:num.grids){
    
    alpha.v  <- alphas[ng, ];
    eta.v <- etas[ng, ];
    
    ln.nu.alpha.C <- D * (lgamma(sum(alpha.v)) - sum(lgamma(alpha.v))); # calculates the constant part of the likilihood ratio
    ln.nu.eta.C <- K * (lgamma(sum(eta.v)) - sum(lgamma(eta.v)));
    
    for (iter in 1:sample.count){
      
      theta.nu.alpha <- sum(t(alpha.v) %*% log(thetas[,,iter])) + ln.nu.alpha.C;
      beta.nu.alpha <- sum(log(betas[,,iter]) %*% eta.v) + ln.nu.eta.C;
      
      nu.alphas[ng, iter] <- theta.nu.alpha + beta.nu.alpha;
      
    }
    
    cat("grid = ", ng, "\n");
    
  } 
  
  nu.alphas;
  
}


calc_likelihood_ratios <- function(log.nu.alphas, base.alpha.idx) {
  
  G <- ncol(log.nu.alphas);
  num.grids <-  nrow(log.nu.alphas);
  ratios <- matrix(0, ncol=1, nrow=num.grids);
  
  
  h1 <- log.nu.alphas[base.alpha.idx, ];    
  for (i in 1:num.grids){
    h <- log.nu.alphas[i, ];
    ratios[i] <- sum(exp(h - h1)) / G;
  }
  
  ratios;   
}