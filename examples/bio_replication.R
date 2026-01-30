GET('https://osf.io/dnc4e/?action=download', write_disk('Study_21_meta_analysis.R', overwrite = TRUE))
source("Study_21_meta_analysis.R")

#### Original experiment

original_A549_climetidine_mu = (A549_control_mean - A549_cimetidine_mean)/
  sqrt(A549_control_sd^2/6 + A549_cimetidine_sd^2/6)

df = (6-1)*(A549_control_sd^2/6 + A549_cimetidine_sd^2/6)^2 / 
  (A549_control_sd^4/6^2 + A549_cimetidine_sd^4/6^2)
  
original_A549_climetidine_pval = 1 - pt(original_A549_climetidine_mu, df = df)


original_ACHN_climetidine_mu = (ACHN_control_mean - ACHN_cimetidine_mean)/
  sqrt(ACHN_control_sd^2/6 + ACHN_cimetidine_sd^2/6)

df = (6-1)*(ACHN_control_sd^2/6 + ACHN_cimetidine_sd^2/6)^2 / 
  (ACHN_control_sd^4/6^2 + ACHN_cimetidine_sd^4/6^2)

original_ACHN_climetidine_pval = 1 - pt(original_ACHN_climetidine_mu, df = df)

original_A549_doxorubicin_mu = (A549_control_mean - A549_dox_mean)/
  sqrt(A549_control_sd^2/6 + A549_dox_sd^2/6)

df = (6-1)*(A549_control_sd^2/6 + A549_dox_sd^2/6)^2 / 
  (A549_control_sd^4/6^2 + A549_dox_sd^4/6^2)

original_A549_doxorubicin_pval = 1 - pt(original_A549_doxorubicin_mu, df = df)


#### Replicated experiment

replicated_A549_climetidine_mu = t.test(A549_C, A549_TxCim)$statistic
replicated_A549_climetidine_pval = t.test(A549_C, A549_TxCim, alternative = "greater")$p.value


replicated_ACHN_climetidine_mu = t.test(ACHN_C, ACHN_TxCim)$statistic
replicated_ACHN_climetidine_pval = t.test(ACHN_C, ACHN_TxCim, alternative = "greater")$p.value

replicated_A549_doxorubicin_mu = t.test(A549_C, A549_TxDox)$statistic
replicated_A549_doxorubicin_pval = t.test(A549_C, A549_TxDox, alternative = "greater")$p.value


#### Analysis

alpha = 0.05

mu1 = c(original_A549_climetidine_mu, original_ACHN_climetidine_mu,
        original_A549_doxorubicin_mu)

pval1 = c(original_A549_climetidine_pval, original_ACHN_climetidine_pval,
          original_A549_doxorubicin_pval)

which(pval1 <= alpha/3)

w = c(0, 0, 1)

pval2 = c(replicated_A549_climetidine_pval, replicated_ACHN_climetidine_pval,
          replicated_A549_doxorubicin_pval)

# Original test decisions
(pval1 <= alpha/3) & (pval2 <= alpha/3)

# New test decisions
(pval1 <= alpha/3) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(3*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)


#### Sensitivity analysis

alpha = 0.05

original_A549_climetidine_mu = 1.1*(A549_control_mean - A549_cimetidine_mean)/
  sqrt(A549_control_sd^2/6 + A549_cimetidine_sd^2/6)

df = (6-1)*(A549_control_sd^2/6 + A549_cimetidine_sd^2/6)^2 / 
  (A549_control_sd^4/6^2 + A549_cimetidine_sd^4/6^2)

original_A549_climetidine_pval = 1 - pt(original_A549_climetidine_mu, df = df)


pval1 = c(original_A549_climetidine_pval, original_ACHN_climetidine_pval,
          original_A549_doxorubicin_pval)

which(pval1 <= alpha/3)


library(nleqslv)


mu.hat = c(original_A549_climetidine_mu, original_A549_doxorubicin_mu)

dslnex = function(w, cstar) {
  y <- rep(0,2)
  y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE))) 
  y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE)))
  y[is.na(y)] = -1/alpha
  
  return(y)
}

xstart = mu.hat/sum(mu.hat)

overall_weight_obj = function(cstar){
  
  fstart = dslnex(xstart, cstar)
  wopt = nleqslv(xstart, dslnex, cstar = cstar, method = 'N',
                 control = list(stepmax = 1, allowSingular = TRUE))$x
  
  return((sum(wopt)-1)^2)
}

copt =  optimise(overall_weight_obj, c(-4*max(mu.hat)-(max(mu.hat)^2)/2,
                                       4*max(mu.hat)-(max(mu.hat)^2)/2))$minimum

opt.result = nleqslv(xstart, dslnex, cstar = copt, method = 'N',
                     control = list(stepmax = 1, allowSingular = TRUE))

w = opt.result$x

w = c(w[1], 0, w[2])

print(w)

pval2 = c(replicated_A549_climetidine_pval, replicated_ACHN_climetidine_pval,
          replicated_A549_doxorubicin_pval)

# Original test decisions
(pval1 <= alpha/3) & (pval2 <= alpha/3)

# New test decisions
(pval1 <= alpha/3) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(3*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)

