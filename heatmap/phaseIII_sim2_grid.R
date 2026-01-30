phase_III_sim = function(treatment.means, alpha = 0.05,
                         nsim = 10^4, w.grid){

  m = 2
  
  PoS = w.PoS = w.matrix = matrix(nrow = nsim, ncol = m)
  # PoS = w.PoS = w.matrix = pow2 = w.pow2 = matrix(nrow = nsim, ncol = m)
  
  mu = treatment.means
  
  for(i in 1:nsim){
    
    T1 = mu.hat = rnorm(m, mean = mu)
    T2 = rnorm(m, mean = mu)
    
    w = rep(0,m)
    
    pval1 = 1 - pnorm(T1)
    
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
    
    pval2 = 1 - pnorm(T2)
    
    PoS[i,] = (pval1 <= alpha/m) & (pval2 <= alpha/m)
    w.PoS[i,] = (pval1 <= alpha/m) & (pval2 <= w*alpha) 
    w.matrix[i,] = w
    
    # pow2[i,] = (pval2 <= alpha/m)
    # w.pow2[i,] = (pval2 <= w*alpha)
    
  }
  
  true.non.nulls = which(treatment.means > 0)

  # Marginal PoS
  mPoS.bonf = colMeans(PoS)
  mPoS.w.bonf = colMeans(w.PoS)
  
  # Disjunctive PoS of trials
  dPoS.bonf = mean(rowSums(PoS[,true.non.nulls, drop = F]) > 0)
  dPoS.w.bonf = mean(rowSums(w.PoS[,true.non.nulls, drop = F]) > 0)
  
  # Mean weights
  w.mean = colMeans(w.matrix)
  
  
  # # Marginal power for trial 2
  # pow2 = colMeans(pow2)
  # w.pow2 = colMeans(w.pow2)
  
  
  # return(list(PoS.bonf = PoS.bonf, PoS.w.bonf = PoS.w.bonf,
  #             dPoS.bonf = dPoS.bonf, dPoS.w.bonf = dPoS.w.bonf,
  #             w.mean = w.mean,
  #             pow2 = pow2, w.pow2 = w.pow2))
  
  return(list(mPoS.bonf = mPoS.bonf, mPoS.w.bonf = mPoS.w.bonf,
              dPoS.bonf = dPoS.bonf, dPoS.w.bonf = dPoS.w.bonf,
              w.mean = w.mean))
  
}