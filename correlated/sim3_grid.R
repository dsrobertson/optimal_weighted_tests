phase_III_sim = function(treatment.means, alpha = 0.05, s=0, w.grid2, w.grid3,
                         rho, nsim = 10^4){
  
  m = length(treatment.means)
  
  PoS = w.PoS = w.matrix = matrix(nrow = nsim, ncol = m)

  mu = treatment.means

  corr.matrix = matrix(rep(rho,3*3), nrow = 3)
  diag(corr.matrix) = 1
  
  
  for(i in 1:nsim){
    
    #T1 = mu.hat = rnorm(m, mean = mu)
    T1 = mu.hat = as.vector(rmvnorm(1, mean = mu, sigma = corr.matrix))
    
    #T2 = rnorm(m, mean = mu*(1+s))
    T2 = as.vector(rmvnorm(1, mean = mu*(1+s), sigma = corr.matrix))
    
    w = rep(0,m)
    
    pval1 = 1 - pnorm(T1)
    
    non.null = which(pval1 <= alpha/m)
    
    mu.hat.select = mu.hat[non.null]
    
    if(length(mu.hat.select) == 3){
      
      eval_f0 <- function( x ) {
        return( (pnorm(qnorm(1-(alpha*x[1])) - mu.hat[1]))* 
                (pnorm(qnorm(1-(alpha*x[2])) - mu.hat[2]))*
                (pnorm(qnorm(1-(alpha*x[3])) - mu.hat[3])) )
      }
      
      obj_val = apply(w.grid3, 1, eval_f0)
      
      w = unlist(w.grid3[which.min(obj_val),])
      
      
    } else if(length(mu.hat.select) == 2){
      
      eval_f0 <- function( x ) {
        return( (pnorm(qnorm(1-(alpha*x[1])) - mu.hat.select[1]))* 
                  (pnorm(qnorm(1-(alpha*x[2])) - mu.hat.select[2])) )
      }
      
      obj_val = apply(w.grid2, 1, eval_f0)
      
      w[non.null] = w.grid2[which.min(obj_val),]
      
    } else if(length(mu.hat.select) == 1) {
      
      w[non.null] = 1

      
    } else {
      w = rep(1/m,m)
    }
    
    w.matrix[i,] = w

    pval2 = 1 - pnorm(T2)
    
    PoS[i,] = (pval1 <= alpha/m) & (pval2 <= alpha/m)
    w.PoS[i,] = (pval1 <= alpha/m) & (pval2 <= w*alpha)
    
  }

  true.non.nulls = which(treatment.means > 0)
  
  # Weights
  w.bonf = colMeans(w.matrix)
  
  # Marginal PoS
  PoS.bonf = colMeans(PoS)
  PoS.w.bonf = colMeans(w.PoS)
  
  # Disjunctive PoS of trials
  dPoS.bonf = mean(rowSums(PoS[,true.non.nulls, drop = F]) > 0)
  dPoS.w.bonf = mean(rowSums(w.PoS[,true.non.nulls, drop = F]) > 0)

  return(list(PoS.bonf = PoS.bonf, PoS.w.bonf = PoS.w.bonf, 
              dPoS.bonf = dPoS.bonf, dPoS.w.bonf = dPoS.w.bonf,
              w.bonf = w.bonf))
    
}