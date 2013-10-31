calc_log_nu_alphas <- function (thetas, alphas) {
  
  sample.count   <- dim(thetas)[3];
  D          <- dim(thetas)[2];
  num.grids  <- ncol(alphas);
  nu.alphas  <- matrix(0, nrow=num.grids, ncol=sample.count);
  
  for (ng in 1:num.grids){
    
    
    alpha.v  <- alphas[, ng];        
    log.nu.alpha.C <- D * (lgamma(sum(alpha.v)) - sum(lgamma(alpha.v))); # calculates the constant part of the likilihood ratio
    
    for (iter in 1:sample.count){
      nu.alphas[ng, iter] <- sum(t(alpha.v) %*% log(thetas[,,iter])) + log.nu.alpha.C;  
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
  
  cat("Computing log nu alphas: \n\t number of grids = ", num.grids);
  cat("\n\t number of samples = ", sample.count);
  cat("\n\t number of topics = ", K);
  cat("\n\t number of docs = ", D);
  cat("\n\t vocabulary size = ", V, "\n\t ");

  log.nu.alphas  <- matrix(0, nrow=num.grids, ncol=sample.count);
  
  for (ng in 1:num.grids){
    
    alpha.v  <- alphas[ng, ]; # gets each grid 
    eta.v <- etas[ng, ]; # gets each grid 
    
    log.nu.alpha.C <- D * (lgamma(sum(alpha.v)) - sum(lgamma(alpha.v))); # calculates the constant part of the likilihood ratio
    log.nu.eta.C <- K * (lgamma(sum(eta.v)) - sum(lgamma(eta.v)));
    
    for (iter in 1:sample.count){
      
      theta.nu.alpha <- sum(t(alpha.v) %*% log(thetas[,,iter])) + log.nu.alpha.C;
      beta.nu.alpha <- sum(log(betas[,,iter]) %*% eta.v) + log.nu.eta.C;
      
      log.nu.alphas[ng, iter] <- theta.nu.alpha + beta.nu.alpha;
      
    }
    if (ng %% 10 == 0) { cat("."); } 
    if (ng %% 1000 == 0) { cat("\n"); } 
    
  } 
  
  cat("\nCompleted.")
  
  log.nu.alphas;
  
}


calc_likelihood_ratios <- function(log.nu.alphas, base.alpha.idx) {
  ## Computes likelihood ratios 
  ##
  ## Arguments: 
  ##  log.nu.alphas - log of nu alphas (and etas)  
  ##  base.alpha.idx - base alpha  (and eta) index 
  ## Returns: 
  ##  ratios - a vector of ratios for all alpha (and eta) grids  
  ##
  
  num.samples <- ncol(log.nu.alphas);
  num.grids <-  nrow(log.nu.alphas);
  ratios <- matrix(0, ncol=1, nrow=num.grids);
  
  
  h1 <- log.nu.alphas[base.alpha.idx, ];    
  for (i in 1:num.grids){
    h <- log.nu.alphas[i, ];
    ratios[i] <- sum(exp(h - h1)) / num.samples;
  }
  
  ratios;   
}