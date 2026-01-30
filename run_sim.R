nsim = 10^5

####### 2 treatments (1) ####### 

source('sim2_grid.R')

#### Grid search
w1 = seq(0,1, by = 0.005)
w2 = 1 - w1
w.grid = cbind(w1, w2)

theta = seq(0.1, 6, length.out = 30)


### Set up parallel computing

library(doParallel)
library(doRNG)
library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)


### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(0,theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid = w.grid)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          SWER.bonf = res$SWER.bonf
                          SWER.w.bonf = res$SWER.w.bonf
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      SWER.bonf = res$SWER.bonf,
                                      SWER.w.bonf = res$SWER.w.bonf
                          ))
                          
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

SWER.bonf = results_comb$SWER.bonf
SWER.w.bonf = results_comb$SWER.w.bonf
SWER = cbind(SWER.bonf, SWER.w.bonf)

save(dPoS, SWER, file = '2_treatments_sim1.Rdata')


#### Plot

library(ggplot2)
library(ggtext)
library(patchwork)

data.dPoS <- data.frame(theta = theta, w.bonf = dPoS[,2], bonf = dPoS[,1])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,6)) +
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1))

# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  theme(legend.position = "none") + ylim(c(0,0.08)) + xlim(c(0,6))

PoS_plots = p1 + p2 + 
  plot_annotation(title = "<b>&theta;</b> = (0, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) +
  plot_layout(nrow = 2)

ggsave('PoS_2T_sim1_revised.pdf', PoS_plots, width = 20, height = 16,
       units = "cm", device = cairo_pdf)



#### SWER under global null

set.seed(7)

res = phase_III_sim(rep(0,2), nsim = 10^6, w.grid = w.grid)

c(res$SWER.bonf, res$SWER.w.bonf)


####### 2 treatments (2) #######

## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)


### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(theta[i]/2, theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid = w.grid)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          mPoS.bonf.H1 = res$PoS.bonf[1]
                          mPoS.w.bonf.H1 = res$PoS.w.bonf[1]
                          
                          mPoS.bonf.H2 = res$PoS.bonf[2]
                          mPoS.w.bonf.H2 = res$PoS.w.bonf[2]
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      mPoS.bonf.H1 = mPoS.bonf.H1,
                                      mPoS.w.bonf.H1 = mPoS.w.bonf.H1,
                                      mPoS.bonf.H2 = mPoS.bonf.H2,
                                      mPoS.w.bonf.H2 = mPoS.w.bonf.H2))
                        }

Sys.time()

stopCluster(cl)


dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

mPoS.bonf.H1 = results_comb$mPoS.bonf.H1
mPoS.w.bonf.H1 = results_comb$mPoS.w.bonf.H1
mPoS.H1 = cbind(mPoS.bonf.H1, mPoS.w.bonf.H1)


mPoS.bonf.H2 = results_comb$mPoS.bonf.H2
mPoS.w.bonf.H2 = results_comb$mPoS.w.bonf.H2
mPoS.H2 = cbind(mPoS.bonf.H2, mPoS.w.bonf.H2)

#save(dPoS, mPoS.H1, mPoS.H2,
#     file = '2_treatments_sim2.Rdata')


#### Disjunctive PoS

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,6)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(str2expression(paste("H[1]"," or ","H[2]", sep = "~" )))

p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.10,0.08)) + xlim(c(0,6))

# PoS_plots = p1 + p2 + plot_layout(nrow = 2)
# 
# ggsave('dPoS_2T_sim2.pdf', PoS_plots, width = 20, height = 16, units = "cm")


#### Marginal PoS

## H1

data.H1 <- data.frame(theta = theta, bonf = mPoS.H1[,1], w.bonf = mPoS.H1[,2])

p1.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("mPoS") + xlab(expression(theta)) + 
  scale_color_manual(values = c("blue", "red")) + xlim(c(0,6)) + 
  labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[1])) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")

# Difference
p2.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  xlab(expression(theta)) + ylab("Difference in mPoS") + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.10,0.08)) + xlim(c(0,6))


## H2

data.H2 <- data.frame(theta = theta, bonf = mPoS.H2[,1], w.bonf = mPoS.H2[,2])

p1.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("mPoS") + xlab(expression(theta)) + xlim(c(0,6)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[2])) + theme(plot.title = element_text(hjust = 0.5))

# Difference
p2.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.10,0.08)) + xlim(c(0,6))


layout_design <- "
  ABC
  DDD
  EFG
"

PoS_plots = p1 + p1.H1 + p1.H2 + 
  guide_area() +
  p2 + p2.H1 + p2.H2 + 
  
  plot_layout(design = layout_design,
              guides = "collect",
              heights = c(1, 0.13, 1)) +
  
  plot_annotation(title = "<b>&theta;</b> = (&theta;/2, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'top', legend.box.margin = margin(b = 20))

ggsave('PoS_2T_sim2.pdf', PoS_plots, width = 28, height = 18,
       units = "cm", device = cairo_pdf)


####### 2 treatments (3) #######

## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)


### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          
                          treatment.means = c(theta[i], theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid = w.grid)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          mPoS.bonf.H1 = res$PoS.bonf[1]
                          mPoS.w.bonf.H1 = res$PoS.w.bonf[1]
                          
                          mPoS.bonf.H2 = res$PoS.bonf[2]
                          mPoS.w.bonf.H2 = res$PoS.w.bonf[2]
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      mPoS.bonf.H1 = mPoS.bonf.H1,
                                      mPoS.w.bonf.H1 = mPoS.w.bonf.H1,
                                      mPoS.bonf.H2 = mPoS.bonf.H2,
                                      mPoS.w.bonf.H2 = mPoS.w.bonf.H2))
                          
                        }

Sys.time()

stopCluster(cl)


dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

mPoS.bonf.H1 = results_comb$mPoS.bonf.H1
mPoS.w.bonf.H1 = results_comb$mPoS.w.bonf.H1
mPoS.H1 = cbind(mPoS.bonf.H1, mPoS.w.bonf.H1)

mPoS.bonf.H2 = results_comb$mPoS.bonf.H2
mPoS.w.bonf.H2 = results_comb$mPoS.w.bonf.H2
mPoS.H2 = cbind(mPoS.bonf.H2, mPoS.w.bonf.H2)

# save(dPoS, mPoS.H1, mPoS.H2,
#      file = '2_treatments_sim3.Rdata')


#### Disjunctive PoS

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,6)) +  
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(str2expression(paste("H[1]"," or ","H[2]", sep = "~" )))


p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.01,0.07)) + xlim(c(0,6))

# PoS_plots = p1 + p2 + plot_layout(nrow = 2)
# 
# ggsave('dPoS_2T_sim3.pdf', PoS_plots, width = 20, height = 16, units = "cm")


#### Marginal PoS

## H1

data.H1 <- data.frame(theta = theta, bonf = mPoS.H1[,1], w.bonf = mPoS.H1[,2])

p1.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("mPoS") + xlab(expression(theta)) + 
  scale_color_manual(values = c("blue", "red")) + xlim(c(0,6)) +
  labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[1])) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")

# Difference
p2.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.01,0.07)) + xlim(c(0,6))


## H2

data.H2 <- data.frame(theta = theta, bonf = mPoS.H2[,1], w.bonf = mPoS.H2[,2])

p1.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("mPoS") + xlab(expression(theta)) + xlim(c(0,6)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[2])) + theme(plot.title = element_text(hjust = 0.5))

# Difference
p2.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(-0.01,0.07)) + xlim(c(0,6))



layout_design <- "
  ABC
  DDD
  EFG
"


PoS_plots = p1 + p1.H1 + p1.H2 + 
  guide_area() +
  p2 + p2.H1 + p2.H2 + 
  
  plot_layout(
    design = layout_design, 
    guides = "collect", 
    # Relative heights: Rows are 1, Legend is 0.1 (thin)
    heights = c(1, 0.13, 1, 1) 
  ) +
  
  plot_annotation(title = "<b>&theta;</b> = (&theta;, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'top', legend.box.margin = margin(b = 20))

ggsave('PoS_2T_sim3.pdf', PoS_plots, width = 28, height = 18,
       units = "cm", device = cairo_pdf)


####### 3 treatments (1) #######
library(ggplot2)
library(cowplot)

# library(nleqslv)

source('sim3_grid.R')

w1 = seq(0,1, by = 0.005)
w2 = 1 - w1
w.grid2 = cbind(w1, w2)

w1 = seq(0,1, by = 0.005)
w2 = seq(0,1, by = 0.005)
w.grid3 = expand.grid(w1, w2)
w.grid3 = w.grid3[which(rowSums(w.grid3) <= 1),]
w.grid3 = cbind(w.grid3, 1 - rowSums(w.grid3))

theta = seq(0.1, 6, length.out = 30)

## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(rep(0,2),theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid2 = w.grid2, w.grid3 = w.grid3)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          SWER.bonf = res$SWER.bonf
                          SWER.w.bonf = res$SWER.w.bonf
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      SWER.bonf = SWER.bonf,
                                      SWER.w.bonf = SWER.w.bonf))
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

SWER.bonf = results_comb$SWER.bonf
SWER.w.bonf = results_comb$SWER.w.bonf
SWER = cbind(SWER.bonf, SWER.w.bonf)

# save(dPoS, SWER, file = '3_treatments_sim1.Rdata')


#### Disjunctive PoS 

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,6)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1))

# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  scale_y_continuous(breaks=seq(0,0.11,0.02), limits = c(0,0.11)) + xlim(c(0,6))

PoS_plots = p1 + p2 + 
  plot_annotation(title = "<b>&theta;</b> = (0, 0, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) +
  plot_layout(nrow = 2)

ggsave('dPoS_3T_sim1.pdf', PoS_plots, width = 20, height = 16,
       units = "cm", device = cairo_pdf)


#### SWER under global null

set.seed(7)

res = phase_III_sim(rep(0,3), nsim = 10^6,
                    w.grid2 = w.grid2, w.grid3 = w.grid3)

c(res$SWER.bonf, res$SWER.w.bonf)


####### 3 treatments (2) #######


## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

theta = seq(0.1, 3, length.out = 15)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)


### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(theta[i]/2,theta[i],2*theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid2 = w.grid2, w.grid3 = w.grid3)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          mPoS.bonf.H1 = res$PoS.bonf[1]
                          mPoS.w.bonf.H1 = res$PoS.w.bonf[1]
                          
                          mPoS.bonf.H2 = res$PoS.bonf[2]
                          mPoS.w.bonf.H2 = res$PoS.w.bonf[2]
                          
                          mPoS.bonf.H3 = res$PoS.bonf[3]
                          mPoS.w.bonf.H3 = res$PoS.w.bonf[3]
                          
                          w.bonf = res$w.bonf
                          
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      mPoS.bonf.H1 = mPoS.bonf.H1,
                                      mPoS.w.bonf.H1 = mPoS.w.bonf.H1,
                                      mPoS.bonf.H2 = mPoS.bonf.H2,
                                      mPoS.w.bonf.H2 = mPoS.w.bonf.H2,
                                      mPoS.bonf.H3 = mPoS.bonf.H3,
                                      mPoS.w.bonf.H3 = mPoS.w.bonf.H3,
                                      w.bonf = w.bonf))
                          
                        }

Sys.time()

stopCluster(cl)

w.bonf = matrix(results_comb$w.bonf, ncol = 3, byrow = TRUE)

dPoS.bonf = matrix(results_comb$dPoS.bonf, ncol = 3, byrow = TRUE)[,1]
dPoS.w.bonf = matrix(results_comb$dPoS.w.bonf, ncol = 3, byrow = TRUE)[,1]
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)


mPoS.bonf.H1 = matrix(results_comb$mPoS.bonf.H1, ncol = 3, byrow = TRUE)[,1]
mPoS.w.bonf.H1 = matrix(results_comb$mPoS.w.bonf.H1, ncol = 3, byrow = TRUE)[,1]
mPoS.H1 = cbind(mPoS.bonf.H1, mPoS.w.bonf.H1)

mPoS.bonf.H2 = matrix(results_comb$mPoS.bonf.H2, ncol = 3, byrow = TRUE)[,1]
mPoS.w.bonf.H2 = matrix(results_comb$mPoS.w.bonf.H2, ncol = 3, byrow = TRUE)[,1]
mPoS.H2 = cbind(mPoS.bonf.H2, mPoS.w.bonf.H2)

mPoS.bonf.H3 = matrix(results_comb$mPoS.bonf.H3, ncol = 3, byrow = TRUE)[,1]
mPoS.w.bonf.H3 = matrix(results_comb$mPoS.w.bonf.H3, ncol = 3, byrow = TRUE)[,1]
mPoS.H3 = cbind(mPoS.bonf.H3, mPoS.w.bonf.H3)


# save(dPoS, mPoS.H1, mPoS.H2, mPoS.H3, w.bonf,
#      file = '3_treatments_sim2.Rdata')


#### Disjunctive PoS 

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,3)) +  
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1))

# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  scale_y_continuous(breaks=seq(0,0.1,0.02), limits = c(0,0.1)) + xlim(c(0,3))

PoS_plots = p1 + p2 + plot_layout(nrow = 2) +
  plot_annotation(title = "<b>&theta;</b> = (&theta;/2, &theta;, 2&theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19)))

ggsave('dPoS_3T_sim2.pdf', PoS_plots, width = 20, height = 16,
       units = "cm", device = cairo_pdf)


#### Marginal PoS

## H1

data.H1 <- data.frame(theta = theta, bonf = mPoS.H1[,1], w.bonf = mPoS.H1[,2],
                      weights = w.bonf[,1])

p1.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("Marginal PoS") + xlab(expression(theta)) + 
  scale_color_manual(values = c("blue", "red")) + xlim(c(0,3)) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[1])) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")

# Difference
p2.H1 <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  #scale_y_continuous(breaks=seq(-0.06,0.10,by=0.02), limits = c(-0.06,0.1)) +
  ylim(c(-0.06, 0.1)) + xlim(c(0,3))

# Weights
p2.H1.weights <- ggplot(data.H1, aes(x = theta)) +
  geom_line(aes(y = weights), color = 'purple') +
  ylab("Mean weight") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(0,1)) + xlim(c(0,3))


## H2

data.H2 <- data.frame(theta = theta, bonf = mPoS.H2[,1], w.bonf = mPoS.H2[,2],
                     weights = w.bonf[,2])

p1.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("") + xlab(expression(theta)) + xlim(c(0,3)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[2])) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")

# Difference
p2.H2 <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  #scale_y_continuous(breaks=seq(-0.06,0.10,by=0.02), limits = c(-0.06,0.1)) +
  ylim(c(-0.06, 0.1)) + xlim(c(0,3))

# Weights
p2.H2.weights <- ggplot(data.H2, aes(x = theta)) +
  geom_line(aes(y = weights), color = 'purple') +
  ylab("") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(0,1)) + xlim(c(0,3))


### H3

data.H3 <- data.frame(theta = theta, bonf = mPoS.H3[,1], w.bonf = mPoS.H3[,2],
                      weights = w.bonf[,3])

p1.H3 <- ggplot(data.H3, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("") + xlab(expression(theta)) + xlim(c(0,3)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[3])) + theme(plot.title = element_text(hjust = 0.5))

# Difference
p2.H3 <- ggplot(data.H3, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  #scale_y_continuous(breaks=seq(-0.06,0.10,by=0.02), limits = c(-0.06,0.1)) + 
  ylim(c(-0.06, 0.1)) + xlim(c(0,3))

# Weights
p2.H3.weights <- ggplot(data.H3, aes(x = theta)) +
  geom_line(aes(y = weights), color = 'purple') +
  ylab("") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  ylim(c(0,1)) + xlim(c(0,3))


# Define the layout visually
# A-C are the first row, D is the legend, E-J are the rest
layout_design <- "
  ABC
  DDD
  EFG
  HIJ
"

PoS_plots = p1.H1 + p1.H2 + p1.H3 +                  # Row 1 (A, B, C)
  guide_area() +                         # Legend Placeholder (D)
  p2.H1 + p2.H2 + p2.H3 +                # Row 2 (E, F, G)
  p2.H1.weights + p2.H2.weights + p2.H3.weights + # Row 3 (H, I, J)
  
  plot_layout(
    design = layout_design, 
    guides = "collect", 
    # Relative heights: Rows are 1, Legend is 0.1 (thin)
    heights = c(1, 0.18, 1, 1) 
  ) +
  
  plot_annotation(
    title = "<b>&theta;</b> = (&theta;/2, &theta;, 2&theta;)",
    theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))
  ) & theme(legend.position = 'top', legend.box.margin = margin(b = 20))


ggsave('mPoS_3T_sim2_revised.pdf', PoS_plots, width = 24, height = 20,
       units = "cm", device = cairo_pdf)


####### 5 treatments (1) #######
library(ggplot2)
library(cowplot)

library(nleqslv)

source('sim5_hybrid.R')

theta = seq(0.1, 6, length.out = 30)

w1 = seq(0,1, by = 0.005)
w2 = 1 - w1

w.grid2 = cbind(w1, w2)

w2 = seq(0,1, by = 0.005)
w.grid = expand.grid(w1, w2)
w.grid = w.grid[which(rowSums(w.grid) <= 1),]
w.grid3 = cbind(w.grid, 1 - rowSums(w.grid))


## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(rep(0,4),theta[i])
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid2 = w.grid2, w.grid3 = w.grid3)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          SWER.bonf = res$SWER.bonf
                          SWER.w.bonf = res$SWER.w.bonf
                          
                          n.NA = res$n.NA
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      SWER.bonf = res$SWER.bonf,
                                      SWER.w.bonf = res$SWER.w.bonf,
                                      n.NA = n.NA))
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

SWER.bonf = results_comb$SWER.bonf
SWER.w.bonf = results_comb$SWER.w.bonf
SWER = cbind(SWER.bonf, SWER.w.bonf)

n.NA = results_comb$n.NA

# save(dPoS, SWER, n.NA, file = '5treatments_sim1.Rdata')

print(n.NA)

#### Disjunctive PoS 

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,6)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1))

# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs( color = NULL) + theme(legend.position = "none") +
  scale_y_continuous(breaks=seq(0,0.14,0.02), limits = c(-0.001,0.14)) +
  xlim(c(0,6))

PoS_plots = p1 + p2 + 
  plot_annotation(title = "<b>&theta;</b> = (0, 0, 0, 0, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) + 
  plot_layout(nrow = 2)

ggsave('PoS_5T_sim1_revised.pdf', PoS_plots, width = 20, height = 16,
       units = "cm", device = cairo_pdf)


#### SWER under global null

set.seed(7)

res = phase_III_sim(rep(0,5), nsim = 10^6,
                    w.grid2 = w.grid2, w.grid3 = w.grid3)

c(res$SWER.bonf, res$SWER.w.bonf)




####### 5 treatments (2) #######

## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

# theta = seq(0.1, 4, length.out = 20)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)
registerDoRNG(7)

### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          treatment.means = c(rep(0,2),rep(theta[i],3))
                          
                          res = phase_III_sim(treatment.means, nsim = nsim,
                                              w.grid2 = w.grid2, w.grid3 = w.grid3)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          SWER.bonf = res$SWER.bonf
                          SWER.w.bonf = res$SWER.w.bonf
                          
                          mPoS.bonf.H5 = res$PoS.bonf[5]
                          mPoS.w.bonf.H5 = res$PoS.w.bonf[5]
                          
                          n.NA = res$n.NA
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      SWER.bonf = res$SWER.bonf,
                                      SWER.w.bonf = res$SWER.w.bonf,
                                      mPoS.bonf.H5 = mPoS.bonf.H5,
                                      mPoS.w.bonf.H5 = mPoS.w.bonf.H5,
                                      n.NA = n.NA))
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

SWER.bonf = results_comb$SWER.bonf
SWER.w.bonf = results_comb$SWER.w.bonf
SWER = cbind(SWER.bonf, SWER.w.bonf)

mPoS.bonf.H5 = results_comb$mPoS.bonf.H5
mPoS.w.bonf.H5 = results_comb$mPoS.w.bonf.H5
mPoS.H5 = cbind(mPoS.bonf.H5, mPoS.w.bonf.H5)

n.NA = results_comb$n.NA

# save(dPoS, SWER, mPoS.H5, n.NA,
#      file = '5treatments_sim2.Rdata')

print(n.NA)

#### Disjunctive PoS 

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,5)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) +
  theme(legend.position = "none") + 
  ggtitle(str2expression(paste("H[3]"," or ","H[4]", " or ", 
                               "H[5]", sep = "~" ))) + 
  theme(plot.title = element_text(hjust = 0.5))

# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs( color = NULL) + theme(legend.position = "none") +
  scale_y_continuous(breaks=seq(0,0.16,0.04), limits = c(-0.005,0.16)) + 
  xlim(c(0,5))

# PoS_plots = p1 + p2 + plot_layout(nrow = 2)
# 
# ggsave('dPoS_5T_sim2.pdf', PoS_plots, width = 20, height = 16, units = "cm")


#### Marginal PoS

### H5

data.mPoS.H5 <- data.frame(theta = theta, bonf = mPoS.H5[,1], w.bonf = mPoS.H5[,2])

p1.H5 <- ggplot(data.mPoS.H5, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("Marginal PoS") + xlab(expression(theta)) + xlim(c(0,5)) +
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[5])) + theme(plot.title = element_text(hjust = 0.5))


# Difference
p2.H5 <- ggplot(data.mPoS.H5, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) +
  labs(color = NULL) + theme(legend.position = "none") + 
  scale_y_continuous(breaks=seq(0,0.16,0.04), limits = c(-0.005,0.16)) +
  xlim(c(0,5))

PoS_plots = p1 + p1.H5 + p2 + p2.H5 +
  plot_annotation(title = "<b>&theta;</b> = (0, 0, &theta;, &theta;, &theta;)",
                  theme = theme(plot.title = element_markdown(hjust = 0.5, size = 19))) +
  plot_layout(nrow = 2)

ggsave('PoS_5T_sim2.pdf', PoS_plots, width = 24, height = 16,
       units = "cm", device = cairo_pdf)



####### 5 treatments (3) #######

source('sim5_hybrid_random.R')

## Set up parallel computation

# library(doParallel)
# library(doRNG)
# library(data.table)

mc = length(theta); cl = makeCluster(mc)

registerDoParallel(cl)

registerDoRNG(7)


### Run simulation

Sys.time()

results_comb = foreach (i = 1:length(theta),
                        .combine = function(x,y)rbindlist(list(x,y)),
                        .packages = 'nleqslv') %dopar% {
                          
                          five.T.mean = theta[i]
                          
                          res = phase_III_sim_5T_random(five.T.mean,
                                                        nsim = nsim,
                                              w.grid2 = w.grid2, 
                                              w.grid3 = w.grid3)
                          
                          dPoS.bonf = res$dPoS.bonf
                          dPoS.w.bonf = res$dPoS.w.bonf
                          
                          mPoS.bonf.H1 = res$PoS.bonf[1]
                          mPoS.w.bonf.H1 = res$PoS.w.bonf[1]
                          
                          mPoS.bonf.H5 = res$PoS.bonf[5]
                          mPoS.w.bonf.H5 = res$PoS.w.bonf[5]
                          
                          n.NA = res$n.NA
                          
                          return(list(dPoS.bonf = dPoS.bonf,
                                      dPoS.w.bonf = dPoS.w.bonf,
                                      mPoS.bonf.H1 = mPoS.bonf.H1,
                                      mPoS.w.bonf.H1 = mPoS.w.bonf.H1,
                                      mPoS.bonf.H5 = mPoS.bonf.H5,
                                      mPoS.w.bonf.H5 = mPoS.w.bonf.H5,
                                      n.NA = n.NA))
                        }

Sys.time()

stopCluster(cl)

dPoS.bonf = results_comb$dPoS.bonf
dPoS.w.bonf = results_comb$dPoS.w.bonf
dPoS = cbind(dPoS.bonf, dPoS.w.bonf)

mPoS.bonf.H1 = results_comb$mPoS.bonf.H1
mPoS.w.bonf.H1 = results_comb$mPoS.w.bonf.H1
mPoS.H1 = cbind(mPoS.bonf.H1, mPoS.w.bonf.H1)

mPoS.bonf.H5 = results_comb$mPoS.bonf.H5
mPoS.w.bonf.H5 = results_comb$mPoS.w.bonf.H5
mPoS.H5 = cbind(mPoS.bonf.H5, mPoS.w.bonf.H5)

n.NA = results_comb$n.NA

save(dPoS, mPoS.H1, mPoS.H5, n.NA,
     file = '5treatments_sim3.Rdata')

print(n.NA)

#### Disjunctive PoS 

data.dPoS <- data.frame(theta = theta, bonf = dPoS[,1], w.bonf = dPoS[,2])

p1 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("dPoS") + xlab(expression(theta)) + xlim(c(0,5)) + 
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  theme(legend.position = "none") + 
  ggtitle(str2expression(paste("Any ", "H[i]", sep = "~" ))) + 
  theme(plot.title = element_text(hjust = 0.5))


# Difference
p2 <- ggplot(data.dPoS, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in dPoS") + xlab(expression(theta)) + 
  labs(color = NULL) + theme(legend.position = "none") + 
  scale_y_continuous(breaks=seq(0,0.16,0.02), limits = c(-0.005,0.16)) +
  xlim(c(0,5))

# PoS_plots = p1 + p2 + plot_layout(nrow = 2)
# 
# ggsave('dPoS_5T_sim3.pdf', PoS_plots, width = 20, height = 16, units = "cm")


#### Marginal PoS

### H5

data.mPoS.H5 <- data.frame(theta = theta, bonf = mPoS.H5[,1], w.bonf = mPoS.H5[,2])

p1.H5 <- ggplot(data.mPoS.H5, aes(x = theta)) +
  geom_line(aes(y = bonf, color = "Bonferroni")) +
  geom_line(aes(y = w.bonf, color = "Weighted Bonferroni")) +
  ylab("Marginal PoS") + xlab(expression(theta)) + xlim(c(0,5)) +
  scale_color_manual(values = c("blue", "red")) + labs(colour = NULL) + 
  guides(colour = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(breaks=seq(0,1,by=0.2), limits = c(0,1)) + 
  ggtitle(expression(H[5])) + theme(plot.title = element_text(hjust = 0.5))


# Difference
p2.H5 <- ggplot(data.mPoS.H5, aes(x = theta)) +
  geom_line(aes(y = (w.bonf - bonf)), color = 'darkgreen') +
  ylab("Difference in mPoS") + xlab(expression(theta)) +
  labs(color = NULL) + theme(legend.position = "none") + 
  scale_y_continuous(breaks=seq(0,0.10,0.02), limits = c(0,0.10)) + 
  xlim(c(0,5))

PoS_plots = p1 + p1.H5 + p2 + p2.H5 + 
  plot_layout(nrow = 2) +
  plot_annotation(title = "<b>&theta;</b> = (&theta;<sub>1</sub>, &theta;<sub>2</sub>, &theta;<sub>3</sub>, &theta;<sub>4</sub>, &theta;<sub>5</sub>) with &theta;<sub>i</sub> &#126; U[0, &theta;]",
                  theme = theme(plot.title = element_markdown(hjust = 0.4, size = 19)))
  
ggsave('PoS_5T_sim3.pdf', PoS_plots, width = 24, height = 16,
       units = "cm", device = cairo_pdf)
