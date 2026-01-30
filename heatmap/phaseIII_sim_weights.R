phase_III_sim_weights = function(treatment.means, alpha = 0.05,
                         nsim = 10^4, w.grid){
  
  m = 2
  
  w.matrix = matrix(nrow = nsim, ncol = m)
  
  mu = treatment.means
  
  for(i in 1:nsim){
    
    mu.hat = rnorm(m, mean = mu)
    
    w = rep(0,m)
    
    pval1 = 1 - pnorm(mu.hat)
    
    non.null = which(pval1 <= alpha/m)
    
    if(length(non.null) == 2){
      
      eval_f0 <- function( x ) {
        return( (pnorm(qnorm(1-(alpha*x[1])) - mu.hat[1]))* 
                  (pnorm(qnorm(1-(alpha*x[2])) - mu.hat[2])) )
      }
      
      obj_val = apply(w.grid, 1, eval_f0)
      
      w = w.grid[which.min(obj_val),]
      
    } else if(length(non.null) == 1) {
      
      w[non.null] = 1
      
    } else {
      w = c(1/m, 1/m)
    }
    
    w.matrix[i,] = w
    
  }
  
  # return(colMeans(w.matrix))
  return(w.matrix)
}