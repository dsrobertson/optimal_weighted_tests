####### 2 treatments heatmap ####### 

source('phaseIII_sim2_grid.R')

###### Run simulation with same treatment means for the two trials ######

#### Grid search
w1 = seq(0,1, by = 0.005)
w2 = 1 - w1
w.grid = cbind(w1, w2)

theta1 = theta2 = seq(0, 6, by = 0.1)
theta = expand.grid(theta1, theta2)
theta = theta[-1,]

### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = 32; cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

Sys.time()

results_comb = foreach (i = 1:dim(theta)[1],
                        .combine = function(x,y)rbindlist(list(x,y))) %dopar% {
                          
                          treatment.means = c(theta[i,1], theta[i,2])
                          
                          res = phase_III_sim(treatment.means, nsim = 10^4,
                                              w.grid = w.grid)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          w.mean = res$w.mean
                          
                          mPoS.bonf = res$mPoS.bonf
                          mPoS.w.bonf = res$mPoS.w.bonf
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      w.mean = w.mean,
                                      mPoS.bonf = mPoS.bonf,
                                      mPoS.w.bonf = mPoS.w.bonf
                          ))
                          
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = matrix(results_comb$dPoS.bonf, ncol = 2, byrow = TRUE)[,1]
dPoS.w.bonf = matrix(results_comb$dPoS.w.bonf, ncol = 2, byrow = TRUE)[,1]
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

w.mean = matrix(results_comb$w.mean, ncol = 2, byrow = TRUE)

mPoS.bonf = matrix(results_comb$mPoS.bonf, ncol = 2, byrow = TRUE)
mPoS.w.bonf = matrix(results_comb$mPoS.w.bonf, ncol = 2, byrow = TRUE)

save(theta, dPoS, w.mean, mPoS.bonf, mPoS.w.bonf,
     file = '2_treatments_heatmap.Rdata')


#### Plot mean weights

library(ggplot2)
library(scico)

data.diff.dPoS <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                             w1 = w.mean[,1])

# Heatmap 
weights_heatmap = ggplot(data.diff.dPoS, aes(theta1, theta2, fill= w1)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_scico(palette = "vik", limits = c(0,1)) + 
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  labs(fill = expression(hat(w)[1])) +
  ggtitle(expression(paste("Mean of empirical weights ", hat(w)[1]))) + 
  theme(plot.title = element_text(hjust = 0.5))


ggsave('2T_weights_revised.pdf', weights_heatmap, width = 14, height = 11, units = "cm")


# Histogram

source('phaseIII_sim_weights.R')

set.seed(7)

w.matrix = phase_III_sim_weights(c(2, 2), w.grid = w.grid, nsim = 10^4)

colMeans(w.matrix)

sum(w.matrix[,1]==0)
sum(w.matrix[,1]==1)


w1_hist = ggplot(data = data.frame(weight = w.matrix[,1]), aes(x = weight)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black") +
  ggtitle(expression(paste(theta, " = (2,2)"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(expression(paste("Empirical weight ", hat(w)[1]))) +
  ylab("Frequency")


ggsave('w1_hist_revised.pdf', w1_hist, width = 15, height = 11, units = "cm")


#### Plot difference in dPoS

library(ggplot2)

data.diff.dPoS <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                             diff.PoS = dPoS[,2] - dPoS[,1])

# Heatmap 
PoS_heatmap = ggplot(data.diff.dPoS, aes(theta1, theta2, fill= diff.PoS)) + 
  geom_tile() + theme_minimal() + scale_fill_viridis_c() + 
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  labs(fill = "dPoS difference") 


ggsave('dPoS_2T_heatmap_revised.pdf', PoS_heatmap, width = 15, height = 11,
       units = "cm")


#### Plot difference in mPoS

data.diff.mPoS.H1 <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                             diff.mPoS.H1 = mPoS.w.bonf[,1] - mPoS.bonf[,1])

data.diff.mPoS.H2 <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                                diff.mPoS.H2 = mPoS.w.bonf[,2] - mPoS.bonf[,2])

# Heatmaps

library(scico)
library(patchwork)

mPoS_heatmap_H1 = ggplot(data.diff.mPoS.H1, 
                     aes(theta1, theta2, fill= diff.mPoS.H1)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[1])) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_scico(palette = "vik", limits = c(-0.13, 0.13)) + 
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  labs(fill = "mPoS difference") 

mPoS_heatmap_H2 = ggplot(data.diff.mPoS.H2, 
                         aes(theta1, theta2, fill= diff.mPoS.H2)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[2])) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_scico(palette = "vik", limits = c(-0.13, 0.13)) + 
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  labs(fill = "mPoS difference") 

mPoS_heatmap = mPoS_heatmap_H1 + mPoS_heatmap_H2 + 
  plot_layout(ncol = 2, nrow = 1) + plot_layout(guides = "collect")

ggsave('mPoS_2T_heatmap_revised.pdf', mPoS_heatmap, width = 23, height = 11,
       units = "cm")



###### Run simulations with different treatment means for the two trials ###### 

source('phaseIII_sim2_grid_robust.R')

#### Grid search
w1 = seq(0,1, by = 0.005)
w2 = 1 - w1
w.grid = cbind(w1, w2)

theta1 = theta2 = seq(0, 6, by = 0.1)
theta = expand.grid(theta1, theta2)
theta = theta[-1,]

### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = 32; cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

Sys.time()

results_comb = foreach (i = 1:dim(theta)[1],
                        .combine = function(x,y)rbindlist(list(x,y))) %dopar% {
                          
                          treatment.means = c(0, theta[i,1])
                          treatment.means2 = c(0, theta[i,2])
                          
                          res_a = phase_III_sim_robust(treatment.means,
                                                       treatment.means2,
                                                       nsim = 10^4,
                                                       w.grid = w.grid)
                          
                          treatment.means = c(theta[i,1]/2, theta[i,1])
                          treatment.means2 = c(theta[i,2]/2, theta[i,2])
                          
                          res_b = phase_III_sim_robust(treatment.means,
                                                       treatment.means2,
                                                       nsim = 10^4,
                                                       w.grid = w.grid)
                          
                          
                          treatment.means = c(theta[i,1], theta[i,1])
                          treatment.means2 = c(theta[i,2], theta[i,2])
                          
                          res_c = phase_III_sim_robust(treatment.means,
                                                       treatment.means2,
                                                       nsim = 10^4,
                                                       w.grid = w.grid)
                          
                          dPoS.bonf_a = res_a$dPoS.bonf
                          dPoS.w.bonf_a = res_a$dPoS.w.bonf
                          PoS.bonf_a = res_a$PoS.bonf
                          PoS.w.bonf_a = res_a$PoS.w.bonf
                          w.mean_a = res_a$w.mean
                          
                          dPoS.bonf_b = res_b$dPoS.bonf
                          dPoS.w.bonf_b = res_b$dPoS.w.bonf
                          PoS.bonf_b = res_b$PoS.bonf
                          PoS.w.bonf_b = res_b$PoS.w.bonf
                          w.mean_b = res_b$w.mean
                          
                          dPoS.bonf_c = res_c$dPoS.bonf
                          dPoS.w.bonf_c = res_c$dPoS.w.bonf
                          PoS.bonf_c = res_c$PoS.bonf
                          PoS.w.bonf_c = res_c$PoS.w.bonf
                          w.mean_c = res_c$w.mean
                          
                          return(list(dPoS.bonf_a = dPoS.bonf_a,
                                      dPoS.w.bonf_a = dPoS.w.bonf_a,
                                      PoS.bonf_a = PoS.bonf_a,
                                      PoS.w.bonf_a = PoS.w.bonf_a,
                                      w.mean_a = w.mean_a,
                                      dPoS.bonf_b = dPoS.bonf_b,
                                      dPoS.w.bonf_b = dPoS.w.bonf_b,
                                      PoS.bonf_b = PoS.bonf_b,
                                      PoS.w.bonf_b = PoS.w.bonf_b,
                                      w.mean_b = w.mean_b,
                                      dPoS.bonf_c = dPoS.bonf_c,
                                      dPoS.w.bonf_c = dPoS.w.bonf_c,
                                      PoS.bonf_c = PoS.bonf_c,
                                      PoS.w.bonf_c = PoS.w.bonf_c,
                                      w.mean_c = w.mean_c
                          ))
                          
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf_a = matrix(results_comb$dPoS.bonf_a, ncol = 2, byrow = TRUE)[,1]
dPoS.w.bonf_a = matrix(results_comb$dPoS.w.bonf_a, ncol = 2, byrow = TRUE)[,1]
dPoS_a = cbind(dPoS.bonf_a, dPoS.w.bonf_a)

mPoS.bonf_a = matrix(results_comb$PoS.bonf_a, ncol = 2, byrow = TRUE)
mPoS.w.bonf_a = matrix(results_comb$PoS.w.bonf_a, ncol = 2, byrow = TRUE)
w.mean_a = matrix(results_comb$w.mean_a, ncol = 2, byrow = TRUE)

dPoS.bonf_b = matrix(results_comb$dPoS.bonf_b, ncol = 2, byrow = TRUE)[,1]
dPoS.w.bonf_b = matrix(results_comb$dPoS.w.bonf_b, ncol = 2, byrow = TRUE)[,1]
dPoS_b = cbind(dPoS.bonf_b, dPoS.w.bonf_b)

mPoS.bonf_b = matrix(results_comb$PoS.bonf_b, ncol = 2, byrow = TRUE)
mPoS.w.bonf_b = matrix(results_comb$PoS.w.bonf_b, ncol = 2, byrow = TRUE)
w.mean_b = matrix(results_comb$w.mean_b, ncol = 2, byrow = TRUE)

dPoS.bonf_c = matrix(results_comb$dPoS.bonf_c, ncol = 2, byrow = TRUE)[,1]
dPoS.w.bonf_c = matrix(results_comb$dPoS.w.bonf_c, ncol = 2, byrow = TRUE)[,1]
dPoS_c = cbind(dPoS.bonf_c, dPoS.w.bonf_c)

mPoS.bonf_c = matrix(results_comb$PoS.bonf_c, ncol = 2, byrow = TRUE)
mPoS.w.bonf_c = matrix(results_comb$PoS.w.bonf_c, ncol = 2, byrow = TRUE)
w.mean_c = matrix(results_comb$w.mean_c, ncol = 2, byrow = TRUE)


save(theta, dPoS_a, mPoS.bonf_a, mPoS.w.bonf_a, w.mean_a,
     dPoS_b, mPoS.bonf_b, mPoS.w.bonf_b, w.mean_b,
     dPoS_c, mPoS.bonf_c, mPoS.w.bonf_c, w.mean_c,
     file = '2_treatments_heatmap_robust.Rdata')


#### Plot difference in dPoS

library(ggplot2)
library(ggtext)

data.diff.dPoS_a <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                               diff.PoS = dPoS_a[,2] - dPoS_a[,1],
                               w1 = w.mean_a[,1])

PoS_heatmap_a = ggplot(data.diff.dPoS_a, aes(theta1, theta2, fill= diff.PoS)) + 
  geom_tile() + theme_minimal() + scale_fill_viridis_c() + 
  xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
  labs(fill = "dPoS difference") + 
  plot_annotation(title = "<b>&theta;</b> = (0, &theta;) and <b>&theta;&apos;</b> = (0, &theta;&apos;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.3, size = 16)))

ggsave('PoS_2T_heatmap_robust_a_revised.pdf', PoS_heatmap_a, width = 15, height = 12,
       units = "cm", device = cairo_pdf)


# weights_heatmap_a = ggplot(data.diff.dPoS_a, aes(theta1, theta2, fill= w1)) + 
#   geom_tile() + theme_minimal() + scale_fill_scico(palette = "vik", limits = c(0,1)) + 
#   xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
#   labs(fill = expression(w[1])) 

data.diff.dPoS_b <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                               diff.PoS = dPoS_b[,2] - dPoS_b[,1],
                               w1 = w.mean_b[,1])

PoS_heatmap_b = ggplot(data.diff.dPoS_b, aes(theta1, theta2, fill= diff.PoS)) + 
  geom_tile() + theme_minimal() + scale_fill_viridis_c() + 
  xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
  labs(fill = "dPoS difference") + 
  plot_annotation(title = "<b>&theta;</b> = (&theta;/2, &theta;) and <b>&theta;&apos;</b> = (&theta;&apos;/2, &theta;&apos;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.3, size = 16)))

ggsave('dPoS_2T_heatmap_robust_b_revised.pdf', PoS_heatmap_b, width = 15, height = 12,
       units = "cm", device = cairo_pdf)


data.diff.dPoS_c <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                               diff.PoS = dPoS_c[,2] - dPoS_c[,1],
                               w1 = w.mean_c[,1])

PoS_heatmap_c = ggplot(data.diff.dPoS_c, aes(theta1, theta2, fill= diff.PoS)) + 
  geom_tile() + theme_minimal() + scale_fill_viridis_c() + 
  xlab(expression(theta)) + ylab(expression(paste(theta, "'"))) +
  labs(fill = "dPoS difference") +
  plot_annotation(title = "<b>&theta;</b> = (&theta;, &theta;) and <b>&theta;&apos;</b> = (&theta;&apos;, &theta;&apos;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.3, size = 16)))

ggsave('dPoS_2T_heatmap_robust_c_revised.pdf', PoS_heatmap_c, width = 15, height = 12,
       units = "cm", device = cairo_pdf)



#### Plot differences in mPoS

library(patchwork)
library(scico)

data.diff.mPoS.H1_b <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                                  diff.mPoS.H1_b = mPoS.w.bonf_b[,1] - mPoS.bonf_b[,1])

data.diff.mPoS.H2_b <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                                  diff.mPoS.H2_b = mPoS.w.bonf_b[,2] - mPoS.bonf_b[,2])

mPoS_heatmap_H1_b = ggplot(data.diff.mPoS.H1_b, 
                           aes(theta1, theta2, fill= diff.mPoS.H1_b)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[1])) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_scico(palette = "vik", limits = c(-0.13, 0.13)) + 
  xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
  labs(fill = "mPoS difference") 

mPoS_heatmap_H2_b = ggplot(data.diff.mPoS.H2_b, 
                           aes(theta1, theta2, fill= diff.mPoS.H2_b)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[2])) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_scico(palette = "vik", limits = c(-0.13, 0.13)) + 
  xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
  labs(fill = "mPoS difference") 

mPoS_heatmap_b = mPoS_heatmap_H1_b + mPoS_heatmap_H2_b + 
  plot_layout(ncol = 2, nrow = 1) + plot_layout(guides = "collect") + 
  plot_annotation(title = "<b>&theta;</b> = (&theta;/2, &theta;) and <b>&theta;&apos;</b> = (&theta;&apos;/2, &theta;&apos;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19)))

ggsave('mPoS_2T_heatmap_robust_b_revised.pdf', mPoS_heatmap_b, width = 23, height = 13,
       units = "cm", device = cairo_pdf)


data.diff.mPoS.H1_c <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                                  diff.mPoS.H1_c = mPoS.w.bonf_c[,1] - mPoS.bonf_c[,1])

data.diff.mPoS.H2_c <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                                  diff.mPoS.H2_c = mPoS.w.bonf_c[,2] - mPoS.bonf_c[,2])

mPoS_heatmap_H1_c = ggplot(data.diff.mPoS.H1_c, 
                           aes(theta1, theta2, fill= diff.mPoS.H1_c)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[1])) + 
  scale_fill_scico(palette = "vik", limits = c(-0.03, 0.03)) + 
  xlab(expression(theta[1])) + ylab(expression(theta[2])) +
  labs(fill = "PoS difference") 

mPoS_heatmap_H2_c = ggplot(data.diff.mPoS.H2_c, 
                           aes(theta1, theta2, fill= diff.mPoS.H2_c)) + 
  geom_tile() + theme_minimal() + ggtitle(expression(H[2])) + 
  scale_fill_scico(palette = "vik", limits = c(-0.03, 0.03)) + 
  xlab(expression(theta[1])) + ylab("") +
  labs(fill = "PoS difference") 

mPoS_heatmap_c = mPoS_heatmap_H1_c + mPoS_heatmap_H2_c + 
  plot_layout(ncol = 2, nrow = 1) + plot_layout(guides = "collect")

ggsave('mPoS_2T_heatmap_robust_c.pdf', mPoS_heatmap_c, width = 23, height = 11,
       units = "cm")


######### Badly misspecified means ########


### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = 32; cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

Sys.time()

results_comb = foreach (i = 1:dim(theta)[1],
                        .combine = function(x,y)rbindlist(list(x,y))) %dopar% {
                          
                          treatment.means = c(theta[i,1]/2, theta[i,1])
                          treatment.means2 = c(theta[i,2], theta[i,2]/2)
                          
                          res = phase_III_sim_robust(treatment.means,
                                                     treatment.means2,
                                                     nsim = 10^4,
                                                     w.grid = w.grid)
                          
                         
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          w.mean = res$w.mean

                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      w.mean = w.mean
                          ))
                          
                        }

Sys.time()

stopCluster(cl)


dPoS.bonf = matrix(results_comb$dPoS.bonf, ncol = 2, byrow = TRUE)[,1]
dPoS.w.bonf = matrix(results_comb$dPoS.w.bonf, ncol = 2, byrow = TRUE)[,1]
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)
w.mean = matrix(results_comb$w.mean, ncol = 2, byrow = TRUE)

save(theta, dPoS, w.mean,file = '2_treatments_heatmap_robust_b.Rdata')


#### Plot difference in dPoS

library(ggplot2)

data.diff.dPoS <- data.frame(theta1 = theta[,1], theta2 = theta[,2],
                               diff.PoS = dPoS[,2] - dPoS[,1],
                               w1 = w.mean[,1])

PoS_heatmap = ggplot(data.diff.dPoS, aes(theta1, theta2, fill= diff.PoS)) + 
  geom_tile() + theme_minimal() + 
  scale_fill_scico(palette = "vik", limits = c(-0.085,0.085)) + 
  xlab(expression(theta)) + ylab(expression(paste(theta,"'"))) +
  labs(fill = "dPoS difference") + 
  plot_annotation(title = "<b>&theta;</b> = (&theta;/2, &theta;) and <b>&theta;&apos;</b> = (&theta;&apos;, &theta;&apos;/2)",
                  theme = theme(plot.title = element_markdown(hjust = 0.3, size = 16)))

ggsave('dPoS_2T_heatmap_robust_revised.pdf', PoS_heatmap, width = 15, height = 12,
       units = "cm", device = cairo_pdf)

###

t1 = 6
t2 = 2

res = phase_III_sim_robust(c(t1/2,t1),
                           c(t2,t2/2),
                           nsim = 10^4,
                           w.grid = w.grid)

res$dPoS.bonf
res$dPoS.w.bonf




