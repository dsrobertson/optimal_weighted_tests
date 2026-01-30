phase_III_sim = function(treatment.means,alpha = 0.05, s=0,
                         nsim = 10^4, w.grid2, w.grid3){
  
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
    
    mu.hat.select = mu.hat[non.null]
    
    if(length(mu.hat.select) == 5){
      
      dslnex = function(w, cstar) {
        y <- rep(0,5)
        
        prob1 = 1 - (alpha * w[1])
        prob2 = 1 - (alpha * w[2])
        prob3 = 1 - (alpha * w[3])
        prob4 = 1 - (alpha * w[4])
        prob5 = 1 - (alpha * w[5])
        
        # Clamp probabilities to be within [0, 1] to avoid NaN in qnorm
        prob1_clamped = max(0, min(1, prob1))
        prob2_clamped = max(0, min(1, prob2))
        prob3_clamped = max(0, min(1, prob3))
        prob4_clamped = max(0, min(1, prob4))
        prob5_clamped = max(0, min(1, prob5))
        
        
        y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob5_clamped) - mu.hat[5], log.p = TRUE))) 
        y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob5_clamped) - mu.hat[5], log.p = TRUE)))
        y[3] <- w[3] - (1/alpha)*pnorm(-mu.hat[3]/2 - (1/mu.hat[3])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob5_clamped) - mu.hat[5], log.p = TRUE)))
        y[4] <- w[4] - (1/alpha)*pnorm(-mu.hat[4]/2 - (1/mu.hat[4])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob5_clamped) - mu.hat[5], log.p = TRUE)))
        y[5] <- w[5] - (1/alpha)*pnorm(-mu.hat[5]/2 - (1/mu.hat[5])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE))) 
        
        return(y)
      }
      
      xstart = mu.hat/sum(mu.hat)
      
      overall_weight_obj = function(cstar){
        
        fstart = dslnex(xstart, cstar)
        wopt = nleqslv(xstart, dslnex, cstar = cstar,
                       control = list())$x
        
        return((sum((wopt))-1)^2 + any(wopt < -1e-10))
      }
      
      copt =  optimise(overall_weight_obj, c(-5*max(mu.hat)-(max(mu.hat)^2)/2,
                                             4+5*max(mu.hat)-(max(mu.hat)^2)/2))$minimum
      
      opt.result = nleqslv(xstart, dslnex, cstar = copt, 
                           control = list())
      
      # if(opt.result$termcd != 1){
      # 
      #   overall_weight_obj2 = function(cstar){
      # 
      #     fstart = dslnex(xstart, cstar)
      #     wopt = nleqslv(xstart, dslnex, cstar = cstar, method = 'N',
      #                    control = list(allowSingular = TRUE))$x
      # 
      #     return((sum(abs(wopt))-1)^2)
      #   }
      # 
      #   copt2 =  optimise(overall_weight_obj2, c(-5*max(mu.hat)-(max(mu.hat)^2)/2,
      #                                            4+5*max(mu.hat)-(max(mu.hat)^2)/2))$minimum
      # 
      #   opt.result2 = nleqslv(xstart, dslnex, cstar = copt2, method = 'N',
      #                         control = list(allowSingular = TRUE))
      # 
      # 
      #   w = opt.result2$x
      # 
      # } else {
      # 
      #   w = opt.result$x
      # }
      
      w = opt.result$x
      
      # if(sum(w) < 0.99 | sum(w) > 1.01){
      #   w = NA
      # }
      
      w = w/sum(w)
      
      if(any(w < -1e-10) | any(w > 1+1e-10)){
        w = NA
      }
      
      
    } else if (length(mu.hat.select) == 4){
      
      mu.hat = mu.hat.select
      
      dslnex = function(w, cstar) {
        y <- rep(0,4)
        
        prob1 = 1 - (alpha * w[1])
        prob2 = 1 - (alpha * w[2])
        prob3 = 1 - (alpha * w[3])
        prob4 = 1 - (alpha * w[4])
        
        # Clamp probabilities to be within [0, 1] to avoid NaN in qnorm
        prob1_clamped = max(0, min(1, prob1))
        prob2_clamped = max(0, min(1, prob2))
        prob3_clamped = max(0, min(1, prob3))
        prob4_clamped = max(0, min(1, prob4))
        
        y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE))) 
        y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE)))
        y[3] <- w[3] - (1/alpha)*pnorm(-mu.hat[3]/2 - (1/mu.hat[3])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob4_clamped) - mu.hat[4], log.p = TRUE)))
        y[4] <- w[4] - (1/alpha)*pnorm(-mu.hat[4]/2 - (1/mu.hat[4])*(cstar - 
                                                                       pnorm(qnorm(prob1_clamped) - mu.hat[1], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob2_clamped) - mu.hat[2], log.p = TRUE) - 
                                                                       pnorm(qnorm(prob3_clamped) - mu.hat[3], log.p = TRUE)))
        
        return(y)
      }
      
      xstart = mu.hat/sum(mu.hat)
      
      overall_weight_obj = function(cstar){
        
        fstart = dslnex(xstart, cstar)
        wopt = nleqslv(xstart, dslnex, cstar = cstar,
                       control = list())$x
        
        return((sum((wopt))-1)^2 + any(wopt < -1e-10))
      }
      
      copt =  optimise(overall_weight_obj, c(-5*max(mu.hat)-(max(mu.hat)^2)/2,
                                             3+5*max(mu.hat)-(max(mu.hat)^2)/2))$minimum
      
      opt.result = nleqslv(xstart, dslnex, cstar = copt, 
                           control = list())
      
      # if(opt.result$termcd != 1){
      # 
      #   overall_weight_obj2 = function(cstar){
      # 
      #     fstart = dslnex(xstart, cstar)
      #     wopt = nleqslv(xstart, dslnex, cstar = cstar, method = 'N',
      #                    control = list(allowSingular = TRUE))$x
      # 
      #     return((sum(wopt)-1)^2)
      #   }
      # 
      #   copt2 =  optimise(overall_weight_obj2, c(-5*max(mu.hat)-(max(mu.hat)^2)/2,
      #                                            3+5*max(mu.hat)-(max(mu.hat)^2)/2))$minimum
      # 
      #   opt.result2 = nleqslv(xstart, dslnex, cstar = copt2, method = 'N',
      #                         control = list(allowSingular = TRUE))
      # 
      #   w.temp = opt.result2$x
      # 
      # } else {
      # 
      #   w.temp = opt.result$x
      # }
      
      w.temp = opt.result$x
      
      w.temp = w.temp/sum(w.temp)
      
      if(any(w < -1e-10) | any(w > 1+1e-10)){
        w = NA
      } else {
        w[non.null] = w.temp
      }
      
    } else if(length(mu.hat.select) == 3){
      
      mu.hat = mu.hat.select
      
      eval_f0 <- function( x ) {
        return( (pnorm(qnorm(1-(alpha*x[1])) - mu.hat[1]))* 
                  (pnorm(qnorm(1-(alpha*x[2])) - mu.hat[2]))*
                  (pnorm(qnorm(1-(alpha*x[3])) - mu.hat[3])) )
      }
      
      obj_val = apply(w.grid3, 1, eval_f0)
      
      w[non.null] = unlist(w.grid3[which.min(obj_val),])
      
    } else if(length(mu.hat.select) == 2){
      
      mu.hat = mu.hat.select
      
      eval_f0 <- function( x ) {
        return( (pnorm(qnorm(1-(alpha*x[1])) - mu.hat[1]))* 
                  (pnorm(qnorm(1-(alpha*x[2])) - mu.hat[2])) )
      }
      
      obj_val = apply(w.grid2, 1, eval_f0)
      
      w[non.null] = unlist(w.grid2[which.min(obj_val),])
      
    } else if(length(mu.hat.select) == 1) {
      
      w[non.null] = 1
      
    } else {
      w = rep(1/m, m)
    }
    
    pval1 = 1 - pnorm(T1)
    pval2 = 1 - pnorm(T2)
    
    PoS[i,] = (pval1 <= alpha/m) & (pval2 <= alpha/m)
    w.PoS[i,] = (pval1 <= alpha/m) & (pval2 <= w*alpha) 
    
  }
  
  # true.nulls = which(treatment.means <= 0)
  true.non.nulls = which(treatment.means > 0)
  
  # Number of phase III trials
  
  # # # Observed means
  # # colMeans(mu.hat.matrix, na.rm = TRUE)
  # # colMeans(mu.hat.select.matrix, na.rm = TRUE)
  # #
  # # # Weights
  # # colMeans(w.matrix, na.rm = TRUE)
  
  n.NA = sum(!complete.cases(w.PoS))
  
  # Marginal PoS
  PoS.bonf = colMeans(PoS, na.rm = TRUE)
  PoS.w.bonf = colMeans(w.PoS, na.rm = TRUE)
  
  # SWER
  true.nulls = which(treatment.means == 0)
  
  if(length(true.nulls) > 0){
    
    SWER.bonf = mean(rowSums(PoS[,true.nulls, drop = F]) > 0)
    SWER.w.bonf = mean(rowSums(w.PoS[,true.nulls, drop = F]) > 0)
    
  } else {SWER.bonf = SWER.w.bonf = rep(NA, nsim)}
  
  # Disjunctive PoS of trials
  dPoS.bonf = mean(rowSums(PoS[,true.non.nulls, drop = F], na.rm = TRUE) > 0)
  dPoS.w.bonf = mean(rowSums(w.PoS[,true.non.nulls, drop = F], na.rm = TRUE) > 0)
  
  return(list(PoS.bonf = PoS.bonf, PoS.w.bonf = PoS.w.bonf,
              dPoS.bonf = dPoS.bonf, dPoS.w.bonf = dPoS.w.bonf,
              SWER.bonf = SWER.bonf, SWER.w.bonf = SWER.w.bonf,
              n.NA = n.NA))
  
}