############################# Phase 3 CAC studies ############################ #

#### Data ####
##### Sample sizes (ITT population with LOCF) #####

# Study 1 (ISTA-BEPO-CSO1)
n1_vehicle = 36
n1_bepreve_10 = 36
n1_bepreve_15 = 35


# Study 2 (CL-S&E-0409071-P)
n2_vehicle = 43
n2_bepreve_10 = 44
n2_bepreve_15 = 43


##### Occular itching #####

### Study 1

## Vehicle vs. Bepreve 1.0%

deg.freedom = n1_vehicle + n1_bepreve_10 - 2

# Visit 3B
occular_pval_study1_bepreve10_3B_3min = 0.0010
occular_pval_study1_bepreve10_3B_5min = 0.0002
occular_pval_study1_bepreve10_3B_7min = 0.0001

occular_mu_study1_bepreve10_3B_3min = qt(1-occular_pval_study1_bepreve10_3B_3min, df = deg.freedom)
occular_mu_study1_bepreve10_3B_5min = qt(1-occular_pval_study1_bepreve10_3B_5min, df = deg.freedom)
occular_mu_study1_bepreve10_3B_7min = qt(1-occular_pval_study1_bepreve10_3B_7min, df = deg.freedom)



# Visit 4
occular_mu_study1_bepreve10_4_3min = 0.9/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve10_4_5min = 1/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve10_4_7min = 1/(0.4/qt(0.975, df = deg.freedom))

occular_pval_study1_bepreve10_4_3min = 1 - pt(occular_mu_study1_bepreve10_4_3min, df = deg.freedom)
occular_pval_study1_bepreve10_4_5min = 1 - pt(occular_mu_study1_bepreve10_4_5min, df = deg.freedom)
occular_pval_study1_bepreve10_4_7min = 1 - pt(occular_mu_study1_bepreve10_4_7min, df = deg.freedom)

# Visit 5
occular_mu_study1_bepreve10_5_3min = 1.3/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve10_5_5min = 1.4/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve10_5_7min = 1.3/(0.4/qt(0.975, df = deg.freedom))

occular_pval_study1_bepreve10_5_3min = 1 - pt(occular_mu_study1_bepreve10_5_3min, df = deg.freedom)
occular_pval_study1_bepreve10_5_5min = 1 - pt(occular_mu_study1_bepreve10_5_5min, df = deg.freedom)
occular_pval_study1_bepreve10_5_7min = 1 - pt(occular_mu_study1_bepreve10_5_7min, df = deg.freedom)


## Vehicle vs. Bepreve 1.5%

deg.freedom = n1_vehicle + n1_bepreve_15 - 2

# Visit 3B
occular_pval_study1_bepreve15_3B_3min = 0.0007
occular_pval_study1_bepreve15_3B_5min = 0.0001
occular_pval_study1_bepreve15_3B_7min = 0.0002

occular_mu_study1_bepreve15_3B_3min = qt(1-occular_pval_study1_bepreve15_3B_3min, df = deg.freedom)
occular_mu_study1_bepreve15_3B_5min = qt(1-occular_pval_study1_bepreve15_3B_5min, df = deg.freedom)
occular_mu_study1_bepreve15_3B_7min = qt(1-occular_pval_study1_bepreve15_3B_7min, df = deg.freedom)


# Visit 4
occular_mu_study1_bepreve15_4_3min = 1.25/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve15_4_5min = 1.45/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve15_4_7min = 1.35/(0.4/qt(0.975, df = deg.freedom))

occular_pval_study1_bepreve15_4_3min = 1 - pt(occular_mu_study1_bepreve15_4_3min, df = deg.freedom)
occular_pval_study1_bepreve15_4_5min = 1 - pt(occular_mu_study1_bepreve15_4_5min, df = deg.freedom)
occular_pval_study1_bepreve15_4_7min = 1 - pt(occular_mu_study1_bepreve15_4_7min , df = deg.freedom)

# Visit 5
occular_mu_study1_bepreve15_5_3min = 1.4/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve15_5_5min = 1.35/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study1_bepreve15_5_7min = 1.25/(0.4/qt(0.975, df = deg.freedom))

occular_pval_study1_bepreve15_5_3min = 1 - pt(occular_mu_study1_bepreve15_5_3min, df = deg.freedom)
occular_pval_study1_bepreve15_5_5min = 1 - pt(occular_mu_study1_bepreve15_5_5min, df = deg.freedom)
occular_pval_study1_bepreve15_5_7min = 1 - pt(occular_mu_study1_bepreve15_5_7min, df = deg.freedom)


### Study 2

## Vehicle vs. Bepreve 1.0%

deg.freedom = n2_vehicle + n2_bepreve_10 - 2

# Visit 3B
occular_pval_study2_bepreve10_3B_3min = 0.0055
occular_pval_study2_bepreve10_3B_5min = 0.0006
occular_pval_study2_bepreve10_3B_7min = 0.0001

occular_mu_study2_bepreve10_3B_3min = qt(1-occular_pval_study2_bepreve10_3B_3min, df = deg.freedom)
occular_mu_study2_bepreve10_3B_5min = qt(1-occular_pval_study2_bepreve10_3B_5min, df = deg.freedom)
occular_mu_study2_bepreve10_3B_7min = qt(1-occular_pval_study2_bepreve10_3B_7min, df = deg.freedom)

# Visit 4
occular_mu_study2_bepreve10_4_3min = 1.2/(0.4/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve10_4_5min = 1.3/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve10_4_7min = 1.2/(0.35/qt(0.975, df = deg.freedom))

occular_pval_study2_bepreve10_4_3min = 1 - pt(occular_mu_study2_bepreve10_4_3min, df = deg.freedom)
occular_pval_study2_bepreve10_4_5min = 1 - pt(occular_mu_study2_bepreve10_4_5min, df = deg.freedom)
occular_pval_study2_bepreve10_4_7min = 1 - pt(occular_mu_study2_bepreve10_4_7min, df = deg.freedom)

# Visit 5
occular_mu_study2_bepreve10_5_3min = 1.4/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve10_5_5min = 1.5/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve10_5_7min = 1.3/(0.4/qt(0.975, df = deg.freedom))

occular_pval_study2_bepreve10_5_3min = 1 - pt(occular_mu_study2_bepreve10_5_3min, df = deg.freedom)
occular_pval_study2_bepreve10_5_5min = 1 - pt(occular_mu_study2_bepreve10_5_5min, df = deg.freedom)
occular_pval_study2_bepreve10_5_7min = 1 - pt(occular_mu_study2_bepreve10_5_7min, df = deg.freedom)


## Vehicle vs. Bepreve 1.5%

deg.freedom = n2_vehicle + n2_bepreve_15 - 2

# Visit 3B
occular_pval_study2_bepreve15_3B_3min = 0.0051
occular_pval_study2_bepreve15_3B_5min = 0.0021
occular_pval_study2_bepreve15_3B_7min = 0.0003

occular_mu_study2_bepreve15_3B_3min = qt(1-occular_pval_study2_bepreve15_3B_3min, df = deg.freedom)
occular_mu_study2_bepreve15_3B_5min = qt(1-occular_pval_study2_bepreve15_3B_5min, df = deg.freedom)
occular_mu_study2_bepreve15_3B_7min = qt(1-occular_pval_study2_bepreve15_3B_7min, df = deg.freedom)


# Visit 4
occular_mu_study2_bepreve15_4_3min = 1.3/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve15_4_5min = 1.3/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve15_4_7min = 1.3/(0.35/qt(0.975, df = deg.freedom))

occular_pval_study2_bepreve15_4_3min = 1 - pt(occular_mu_study2_bepreve15_4_3min, df = deg.freedom)
occular_pval_study2_bepreve15_4_5min = 1 - pt(occular_mu_study2_bepreve15_4_5min, df = deg.freedom)
occular_pval_study2_bepreve15_4_7min = 1 - pt(occular_mu_study2_bepreve15_4_7min, df = deg.freedom)


# Visit 5
occular_mu_study2_bepreve15_5_3min = 1.5/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve15_5_5min = 1.6/(0.35/qt(0.975, df = deg.freedom))
occular_mu_study2_bepreve15_5_7min = 1.4/(0.35/qt(0.975, df = deg.freedom))

occular_pval_study2_bepreve15_5_3min = 1 - pt(occular_mu_study2_bepreve15_5_3min, df = deg.freedom)
occular_pval_study2_bepreve15_5_5min = 1 - pt(occular_mu_study2_bepreve15_5_5min, df = deg.freedom)
occular_pval_study2_bepreve15_5_7min = 1 - pt(occular_mu_study2_bepreve15_5_7min, df = deg.freedom)


##### Conjunctival Redness #####

### Study 1


## Vehicle vs. Bepreve 1.0%

deg.freedom = n1_vehicle + n1_bepreve_10 - 2

# Visit 3B
redness_pval_study1_bepreve10_3B_7min = 0.0148
redness_pval_study1_bepreve10_3B_15min = 0.0549
redness_pval_study1_bepreve10_3B_20min = 0.1093

redness_mu_study1_bepreve10_3B_7min = qt(1 - redness_pval_study1_bepreve10_3B_7min, df = deg.freedom)
redness_mu_study1_bepreve10_3B_15min = qt(1 - redness_pval_study1_bepreve10_3B_15min, df = deg.freedom)
redness_mu_study1_bepreve10_3B_20min = qt(1 - redness_pval_study1_bepreve10_3B_20min, df = deg.freedom)

# Visit 4
redness_pval_study1_bepreve10_4_7min = 0.0074
redness_pval_study1_bepreve10_4_15min = 0.0335
redness_pval_study1_bepreve10_4_20min = 0.0544

redness_mu_study1_bepreve10_4_7min = qt(1 - redness_pval_study1_bepreve10_4_7min, df = deg.freedom)
redness_mu_study1_bepreve10_4_15min = qt(1 - redness_pval_study1_bepreve10_4_15min, df = deg.freedom)
redness_mu_study1_bepreve10_4_20min = qt(1 - redness_pval_study1_bepreve10_4_20min, df = deg.freedom)

# Visit 5
redness_pval_study1_bepreve10_5_7min = 0.00004 # only know <0.0001
redness_pval_study1_bepreve10_5_15min = 0.0001
redness_pval_study1_bepreve10_5_20min = 0.0009


redness_mu_study1_bepreve10_5_7min = qt(1 - redness_pval_study1_bepreve10_5_7min, df = deg.freedom)
redness_mu_study1_bepreve10_5_15min = qt(1 - redness_pval_study1_bepreve10_5_15min, df = deg.freedom)
redness_mu_study1_bepreve10_5_20min = qt(1 - redness_pval_study1_bepreve10_5_20min, df = deg.freedom)



## Vehicle vs. Bepreve 1.5%

deg.freedom = n1_vehicle + n1_bepreve_15 - 2

# Visit 3B
redness_pval_study1_bepreve15_3B_7min = 0.3564
redness_pval_study1_bepreve15_3B_15min = 0.4171
redness_pval_study1_bepreve15_3B_20min = 0.7422

redness_mu_study1_bepreve15_3B_7min = qt(1 - redness_pval_study1_bepreve15_3B_7min, df = deg.freedom)
redness_mu_study1_bepreve15_3B_15min = qt(1 - redness_pval_study1_bepreve15_3B_15min, df = deg.freedom)
redness_mu_study1_bepreve15_3B_20min = qt(1 - redness_pval_study1_bepreve15_3B_20min, df = deg.freedom)

# Visit 4
redness_pval_study1_bepreve15_4_7min = 0.0353
redness_pval_study1_bepreve15_4_15min = 0.0626
redness_pval_study1_bepreve15_4_20min = 0.0953

redness_mu_study1_bepreve15_4_7min = qt(1 - redness_pval_study1_bepreve15_4_7min, df = deg.freedom)
redness_mu_study1_bepreve15_4_15min = qt(1 - redness_pval_study1_bepreve15_4_15min, df = deg.freedom)
redness_mu_study1_bepreve15_4_20min = qt(1 - redness_pval_study1_bepreve15_4_20min, df = deg.freedom)

# Visit 5
redness_pval_study1_bepreve15_5_7min = 0.0011
redness_pval_study1_bepreve15_5_15min = 0.0061
redness_pval_study1_bepreve15_5_20min = 0.0482

redness_mu_study1_bepreve15_5_7min = qt(1 - redness_pval_study1_bepreve15_5_7min, df = deg.freedom)
redness_mu_study1_bepreve15_5_15min = qt(1 - redness_pval_study1_bepreve15_5_15min, df = deg.freedom)
redness_mu_study1_bepreve15_5_20min = qt(1 - redness_pval_study1_bepreve15_5_20min, df = deg.freedom)


### Study 2

## Vehicle vs. Bepreve 1.0%

deg.freedom = n2_vehicle + n2_bepreve_10 - 2

# Visit 3B
redness_pval_study2_bepreve10_3B_7min = 0.0053
redness_pval_study2_bepreve10_3B_15min = 0.0168
redness_pval_study2_bepreve10_3B_20min = 0.0407

redness_mu_study2_bepreve10_3B_7min = qt(1 - redness_pval_study2_bepreve10_3B_7min, df = deg.freedom)
redness_mu_study2_bepreve10_3B_15min = qt(1 - redness_pval_study2_bepreve10_3B_15min, df = deg.freedom)
redness_mu_study2_bepreve10_3B_20min = qt(1 - redness_pval_study2_bepreve10_3B_20min, df = deg.freedom)


# Visit 4
redness_pval_study2_bepreve10_4_7min = 0.0006
redness_pval_study2_bepreve10_4_15min = 0.0356
redness_pval_study2_bepreve10_4_20min = 0.1026

redness_mu_study2_bepreve10_4_7min = qt(1 - redness_pval_study2_bepreve10_4_7min, df = deg.freedom)
redness_mu_study2_bepreve10_4_15min = qt(1 - redness_pval_study2_bepreve10_4_15min, df = deg.freedom)
redness_mu_study2_bepreve10_4_20min = qt(1 - redness_pval_study2_bepreve10_4_20min, df = deg.freedom)

# Visit 5
redness_pval_study2_bepreve10_5_7min = 0.0001
redness_pval_study2_bepreve10_5_15min = 0.0020
redness_pval_study2_bepreve10_5_20min = 0.1485

redness_mu_study2_bepreve10_5_7min = qt(1 - redness_pval_study2_bepreve10_5_7min, df = deg.freedom)
redness_mu_study2_bepreve10_5_15min = qt(1 - redness_pval_study2_bepreve10_5_15min, df = deg.freedom)
redness_mu_study2_bepreve10_5_20min = qt(1 - redness_pval_study2_bepreve10_5_20min, df = deg.freedom)


## Vehicle vs. Bepreve 1.5%

deg.freedom = n2_vehicle + n2_bepreve_15 - 2

# Visit 3B
redness_pval_study2_bepreve15_3B_7min = 0.5472
redness_pval_study2_bepreve15_3B_15min = 0.3882
redness_pval_study2_bepreve15_3B_20min = 0.5

redness_mu_study2_bepreve15_3B_7min = qt(1 - redness_pval_study2_bepreve15_3B_7min, df = deg.freedom)
redness_mu_study2_bepreve15_3B_15min = qt(1 - redness_pval_study2_bepreve15_3B_15min, df = deg.freedom)
redness_mu_study2_bepreve15_3B_20min = qt(1 - redness_pval_study2_bepreve15_3B_20min, df = deg.freedom)

# Visit 4
redness_pval_study2_bepreve15_4_7min = 0.1067
redness_pval_study2_bepreve15_4_15min = 0.3598
redness_pval_study2_bepreve15_4_20min = 0.5909

redness_mu_study2_bepreve15_4_7min = qt(1 - redness_pval_study2_bepreve15_4_7min, df = deg.freedom)
redness_mu_study2_bepreve15_4_15min = qt(1 - redness_pval_study2_bepreve15_4_15min, df = deg.freedom)
redness_mu_study2_bepreve15_4_20min = qt(1 - redness_pval_study2_bepreve15_4_20min, df = deg.freedom)

# Visit 5
redness_pval_study2_bepreve15_5_7min = 0.0031
redness_pval_study2_bepreve15_5_15min = 0.0114
redness_pval_study2_bepreve15_5_20min = 0.2251

redness_mu_study2_bepreve15_5_7min = qt(1 - redness_pval_study2_bepreve15_5_7min, df = deg.freedom)
redness_mu_study2_bepreve15_5_15min = qt(1 - redness_pval_study2_bepreve15_5_15min, df = deg.freedom)
redness_mu_study2_bepreve15_5_20min = qt(1 - redness_pval_study2_bepreve15_5_20min, df = deg.freedom)


#### Analysis ####


##### Visit 3B, 7 min post-CAC #####

m = 4
alpha = 0.05


# Weights based on study 

pval1 = c(occular_pval_study1_bepreve10_3B_7min, occular_pval_study1_bepreve15_3B_7min,
          redness_pval_study1_bepreve10_3B_7min, redness_pval_study1_bepreve15_3B_7min)

which(pval1 < alpha/m)


mu.hat = c(occular_mu_study1_bepreve10_3B_7min, occular_mu_study1_bepreve15_3B_7min)

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
# sum(w)

w = c(w, 0, 0)


pval2 = c(occular_pval_study2_bepreve10_3B_7min, occular_pval_study2_bepreve15_3B_7min,
          redness_pval_study2_bepreve10_3B_7min, redness_pval_study2_bepreve15_3B_7min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)



##### Visit 4, 7 min post-CAC #####

m = 4
alpha = 0.05


# Weights based on study 

pval1 = c(occular_pval_study1_bepreve10_4_7min, occular_pval_study1_bepreve15_4_7min,
          redness_pval_study1_bepreve10_4_7min, redness_pval_study1_bepreve15_4_7min)

which(pval1 < alpha/m)

mu.hat = c(occular_mu_study1_bepreve10_4_7min, occular_mu_study1_bepreve15_4_7min,
           redness_mu_study1_bepreve10_4_7min)

dslnex = function(w, cstar) {
  y <- rep(0,3)
  y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE))) 
  y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE)))
  y[3] <- w[3] - (1/alpha)*pnorm(-mu.hat[3]/2 - (1/mu.hat[3])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE)))
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
# sum(w)

w = c(w, 0)


pval2 = c(occular_pval_study2_bepreve10_4_7min, occular_pval_study2_bepreve15_4_7min,
          redness_pval_study2_bepreve10_4_7min, redness_pval_study2_bepreve15_4_7min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)



##### Visit 5, 7 min post-CAC #####

library(nleqslv)

m = 4
alpha = 0.05


# Weights based on study 1
pval1 = c(occular_pval_study1_bepreve10_5_7min, occular_pval_study1_bepreve15_5_7min,
          redness_pval_study1_bepreve10_5_7min, redness_pval_study1_bepreve15_5_7min)

which(pval1 < alpha/m )

mu.hat = c(occular_mu_study1_bepreve10_5_7min, occular_mu_study1_bepreve15_5_7min,
           redness_mu_study1_bepreve10_5_7min, redness_mu_study1_bepreve15_5_7min)

dslnex = function(w, cstar) {
  y <- rep(0,4)
  y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE))) 
  y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE)))
  y[3] <- w[3] - (1/alpha)*pnorm(-mu.hat[3]/2 - (1/mu.hat[3])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE)))
  y[4] <- w[4] - (1/alpha)*pnorm(-mu.hat[4]/2 - (1/mu.hat[4])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE)))
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
# sum(w)


pval2 = c(occular_pval_study2_bepreve10_5_7min, occular_pval_study2_bepreve15_5_7min,
          redness_pval_study2_bepreve10_5_7min, redness_pval_study2_bepreve15_5_7min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= w*alpha)


# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)




##### Conjunctival Redness, Visit 3B #####

m = 6
alpha = 0.05

# Weights based on study 


pval1 = c(redness_pval_study1_bepreve10_3B_7min, redness_pval_study1_bepreve10_3B_15min,
          redness_pval_study1_bepreve10_3B_20min,
          redness_pval_study1_bepreve15_3B_7min, redness_pval_study1_bepreve15_3B_15min,
          redness_pval_study1_bepreve15_3B_20min)

which(pval1 < alpha/m)


pval2 = c(redness_pval_study2_bepreve10_3B_7min, redness_pval_study2_bepreve10_3B_15min,
          redness_pval_study2_bepreve10_3B_20min,
          redness_pval_study2_bepreve15_3B_7min, redness_pval_study2_bepreve15_3B_15min,
          redness_pval_study2_bepreve15_3B_20min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(m*pval2,1)


#### Hypothetical example of decisions changing

alpha = 0.025

pval1 = c(0.004, redness_pval_study1_bepreve10_3B_15min,
          redness_pval_study1_bepreve10_3B_20min,
          redness_pval_study1_bepreve15_3B_7min, redness_pval_study1_bepreve15_3B_15min,
          redness_pval_study1_bepreve15_3B_20min)

which(pval1 < alpha/m)


w = c(1, rep(0,5))

pval2 = c(redness_pval_study2_bepreve10_3B_7min, redness_pval_study2_bepreve10_3B_15min,
          redness_pval_study2_bepreve10_3B_20min,
          redness_pval_study2_bepreve15_3B_7min, redness_pval_study2_bepreve15_3B_15min,
          redness_pval_study2_bepreve15_3B_20min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)


##### Conjunctival Redness, Visit 4 #####

m = 6
alpha = 0.05


# Weights based on study 

pval1 = c(redness_pval_study1_bepreve10_4_7min, redness_pval_study1_bepreve10_4_15min,
          redness_pval_study1_bepreve10_4_20min,
          redness_pval_study1_bepreve15_4_7min, redness_pval_study1_bepreve15_4_15min,
          redness_pval_study1_bepreve15_4_20min)

which(pval1 < alpha/m)

w = c(1, rep(0,5))

pval2 = c(redness_pval_study2_bepreve10_4_7min, redness_pval_study2_bepreve10_4_15min,
          redness_pval_study2_bepreve10_4_20min,
          redness_pval_study2_bepreve15_4_7min, redness_pval_study2_bepreve15_4_15min,
          redness_pval_study2_bepreve15_4_20min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)


##### Conjunctival Redness, Visit 5 #####

m = 6
alpha = 0.05


# Weights based on study 

pval1 = c(redness_pval_study1_bepreve10_5_7min, redness_pval_study1_bepreve10_5_15min,
          redness_pval_study1_bepreve10_5_20min,
          redness_pval_study1_bepreve15_5_7min, redness_pval_study1_bepreve15_5_15min,
          redness_pval_study1_bepreve15_5_20min)


which(pval1 < alpha/m)


mu.hat = c(redness_mu_study1_bepreve10_5_7min, redness_mu_study1_bepreve10_5_15min,
           redness_mu_study1_bepreve10_5_20min,
           redness_mu_study1_bepreve15_5_7min, redness_mu_study1_bepreve15_5_15min)


dslnex = function(w, cstar) {
  y <- rep(0,5)
  y[1] <- w[1] - (1/alpha)*pnorm(-mu.hat[1]/2 - (1/mu.hat[1])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[5])) - mu.hat[5], log.p = TRUE))) 
  y[2] <- w[2] - (1/alpha)*pnorm(-mu.hat[2]/2 - (1/mu.hat[2])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[5])) - mu.hat[5], log.p = TRUE)))
  y[3] <- w[3] - (1/alpha)*pnorm(-mu.hat[3]/2 - (1/mu.hat[3])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[5])) - mu.hat[5], log.p = TRUE)))
  y[4] <- w[4] - (1/alpha)*pnorm(-mu.hat[4]/2 - (1/mu.hat[4])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[5])) - mu.hat[5], log.p = TRUE)))
  y[5] <- w[5] - (1/alpha)*pnorm(-mu.hat[5]/2 - (1/mu.hat[5])*(cstar - 
                                                                 pnorm(qnorm(1-(alpha*w[1])) - mu.hat[1], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[2])) - mu.hat[2], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[3])) - mu.hat[3], log.p = TRUE) - 
                                                                 pnorm(qnorm(1-(alpha*w[4])) - mu.hat[4], log.p = TRUE)))
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
# sum(w)

w = c(w, 0)

pval2 = c(redness_pval_study2_bepreve10_5_7min, redness_pval_study2_bepreve10_5_15min,
          redness_pval_study2_bepreve10_5_20min,
          redness_pval_study2_bepreve15_5_7min, redness_pval_study2_bepreve15_5_15min,
          redness_pval_study2_bepreve15_5_20min)

# Original test decisions
(pval1 <= alpha/m) & (pval2 <= alpha/m)

# New test decisions
(pval1 <= alpha/m) & (pval2 <= w*alpha)

# Original adjusted p-values for study 2
pmin(m*pval2,1)

# New adjusted p-values for study 2
pmin(pval2/w,1)



