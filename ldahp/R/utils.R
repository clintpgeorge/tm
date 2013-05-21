################################################################################################
## Utility Functions
################################################################################################

calc_topic_counts <- function(Z, K)
{
  ## Calculates topic counts 
  # 
  # Inputs: 
  #  Z		- topic selection for each word instance
  #	K		- total topics in a corpus 
  # 
  
  Nt <- array(0, c(1, K)); 
  for (k in 1:K) Nt[k] <- sum(Z == k);
  
  return(Nt); 	
}


normalize <- function (XX, dim=1)
{
  ## Normalizes a given matrix XX 
  # 
  # Inputs: 
  #	XX	- data 
  #	dim - normalizing dimension [1] - column wise, [2] - row wise    
  
  
  if (dim == 1){
    cs <- colSums(XX);
    on <- array(1, c(1, dim(XX)[1]));
    Y <- XX / t(cs %*% on);
  }
  else {
    rs <- rowSums(XX);
    on <- array(1, c(1, dim(XX)[2]));
    Y <- XX / (rs %*% on);
  }
  
  return(Y);
  
}


gen_meshgrid <- function(start, end, interval=0.2){    
  
  # Generates the mesh grid 
  
  meshgrid <- function(a,b) {
    list(
      x=outer(b*0,a,FUN="+"),
      y=outer(b,a*0,FUN="+")
    )
  } 
  
  l <- meshgrid(seq(start, end, by=interval), seq(start, end, by=interval));
  
  alphas <- rbind(as.vector(l$x), as.vector(l$y));
  
}


display_ratios <- function(ratios, start, end, interval, xlabel="alpha-1", ylabel="alpha-2", zlabel="lilelihood ratios"){
  
  # displays the likelihood ratios in a grid 
  
  x <- seq(start, end, by=interval); # e.g. (0.4, 12, 0.4)
  y <- x;
  L <- length(x);
  
  z <- array(0, dim=c(L, L));
  count <- 1;
  for (i in 1:L){
    for (j in 1:L){
      z[i,j] <- ratios[count]; 
      count <- count + 1;
    }
  }
  
  op <- par(bg = "white")
  persp(x, y, z, col = "lightblue")
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
        ltheta = 120, shade = 0.75, ticktype = "detailed",
        xlab = xlabel, ylab = ylabel, zlab = zlabel 
  ) -> res
  
  round(res, 3);
  
}


################################################################################################


compute_thetas <- function(did, Z, K, D, base.alpha.v)
{
  # Re-computes theta for each Gibbs iteration from 
  # the stored z values. This is useful when you have 
  # a dataset with lots of documents and you have 
  # minimal memory 
  
  total.N           <- dim(Z)[1];
  sample.count      <- dim(Z)[2];
  thetas            <- array(0, dim=c(K, D, sample.count));
  
  for (iter in 1:sample.count) {
    zid             <- Z[, iter];
    # base.alpha.v is used because in the Gibbs sampler we use this 
    theta           <- kronecker(matrix(1, 1, D), base.alpha.v); 
    for (i in 1:total.N) { theta[zid[i], did[i]] <- theta[zid[i], did[i]] + 1; }
    for (d in 1:D) { theta[, d] <- rdirichlet(1, theta[, d]); }
    thetas[,,iter]  <- theta;
  }
  
  thetas # return object 
}



compute_thetas_betas <- function(did, wid, Z, K, D, V, base.alpha.v, base.eta)
{
  # Re-computes theta and beta for each Gibbs iteration from 
  # the stored z values. This is useful when you have 
  # a dataset with lots of documents and you have 
  # minimal memory 
  
  total.N         <- dim(Z)[1];
  sample.count    <- dim(Z)[2];
  thetas          <- array(0, dim=c(K, D, sample.count));
  betas           <- array(0, dim=c(K, V, sample.count));
  
  for (iter in 1:sample.count) {
    zid             <- Z[, iter];
    
    # base.alpha.v and base.eta is used, because the Gibbs sampler uses the same 
    theta       <- kronecker(matrix(1, 1, D), base.alpha.v); 
    beta        <- matrix(base.eta, nrow=K, ncol=V);
    
    for (i in 1:total.N) { beta[zid[i], wid[i]] <- beta[zid[i], wid[i]] + 1; }
    for (k in 1:K) { beta[k,] <- rdirichlet(1, beta[k,]); } 
    betas[,,iter]  <- beta;
    
    for (i in 1:total.N) { theta[zid[i], did[i]] <- theta[zid[i], did[i]] + 1; }
    for (d in 1:D) { theta[, d] <- rdirichlet(1, theta[, d]); }
    thetas[,,iter]  <- theta;
    
  }
  
  list(thetas=thetas, betas=betas)
}



