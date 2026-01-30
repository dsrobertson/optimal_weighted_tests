phase_III_sim_compare = function(treatment.means, alpha = 0.05, s=0,
                         nsim = 10^3, w.grid){
  
  m = length(treatment.means)
  
  PoS = w.PoS = h.PoS = matrix(nrow = nsim, ncol = m)
  #mu.hat.matrix = pval.matrix = matrix(nrow = nsim, ncol = m)
  
  mu = treatment.means

  for(i in 1:nsim){
    
    T1 = mu.hat = rnorm(m, mean = mu)
    T2 = rnorm(m, mean = mu*(1+s))
    
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
    
  }
  
  # true.nulls = which(treatment.means <= 0)
  true.non.nulls = which(treatment.means > 0)
  
  R1_11 = sum(PoS[,1] & w.PoS[,1])
  R1_10 = sum(PoS[,1] & !(w.PoS[,1]))
  R1_01 = sum(!(PoS[,1]) & w.PoS[,1])
  R1_00 = sum(!(PoS[,1]) & !(w.PoS[,1]))
  
  R2_11 = sum(PoS[,2] & w.PoS[,2])
  R2_10 = sum(PoS[,2] & !(w.PoS[,2]))
  R2_01 = sum(!(PoS[,2]) & w.PoS[,2])
  R2_00 = sum(!(PoS[,2]) & !(w.PoS[,2]))
  
  R1 = c(R1_11, R1_10, R1_01, R1_00)
  R2 = c(R2_11, R2_10, R2_01, R2_00)
  
  # Marginal PoS
  PoS.bonf = colMeans(PoS)
  PoS.w.bonf = colMeans(w.PoS)
  
  # Disjunctive PoS of trials
  dPoS.bonf = mean(rowSums(PoS[,true.non.nulls, drop = F]) > 0)
  dPoS.w.bonf = mean(rowSums(w.PoS[,true.non.nulls, drop = F]) > 0)
  
  return(list(R1 = R1, R2 = R2, PoS.bonf = PoS.bonf, PoS.w.bonf = PoS.w.bonf,
              dPoS.bonf = dPoS.bonf, dPoS.w.bonf = dPoS.w.bonf))
  
}