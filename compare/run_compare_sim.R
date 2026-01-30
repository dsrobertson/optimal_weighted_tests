source('sim2_compare.R')

#### Grid search
w1 = seq(0,2, by = 0.005)
w2 = 2 - w1
w.grid = cbind(w1, w2)


### Scenario 1

set.seed(7)

theta = 3
theta = c(0, theta)


res = phase_III_sim_compare(theta, alpha = 0.05, s=0,
                            nsim = 10^3, w.grid)

matrix(res$R1, nrow = 2, byrow = TRUE)
matrix(res$R2, nrow = 2, byrow = TRUE)


### Scenario 2

set.seed(7)

theta = 3
theta = c(theta/2, theta)


res = phase_III_sim_compare(theta, alpha = 0.05, s=0,
                            nsim = 10^3, w.grid)

matrix(res$R1, nrow = 2, byrow = TRUE)
matrix(res$R2, nrow = 2, byrow = TRUE)


### Scenario 3

set.seed(7)

theta = 3
theta = c(theta, theta)


res = phase_III_sim_compare(theta, alpha = 0.05, s=0,
                            nsim = 10^3, w.grid)

matrix(res$R1, nrow = 2, byrow = TRUE)
matrix(res$R2, nrow = 2, byrow = TRUE)

