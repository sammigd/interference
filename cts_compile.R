library(here)
library(tidyverse)

n_clus = 200
n_i <- 15
n_time <- 2
n_gamma <- 11
concor = 0.5
n_rep = 1

ngam = 11
gl = seq(from = -.8, to = .8, length.out = 10)
gamma_list = rbind(rep(0, 10),
                   c(gl))
gamma_cts = cbind(gamma_list, rep(0, 2))


here()

alliters = list.files(path = here('ch3truth'))

truth = array(NA, dim = c(n_clus, n_time, n_gamma, length(alliters))) ##clus, time, g


for (rep in 1:length(alliters)){
  load(paste0('ch3truth/', alliters[rep]))
  truth[,,,rep] <- clust_avg_y[,,,1]
}

avg_truth = apply(truth, c(2,3), mean, na.rm = T)
#############################################
alliters_sim = list.files(path = 'ch3sim')
alliters_wts = list.files(path = 'ch3sim_wts')
alliters_sum = list.files(path = 'ch3sim_num')


y_haj = array(NA, dim = c(200, n_time, n_gamma, length(alliters))) ##clus, time, g
wt_haj =  array(NA, dim = c(200, n_time, n_gamma, length(alliters))) 
num =  array(NA, dim = c(200, n_time, n_gamma, length(alliters))) 


for (rep in 1:length(alliters_sim)){
  load(paste0('ch3sim/', alliters_sim[rep]))
  load(paste0('ch3sim_wts/', alliters_wts[rep]))
  load(paste0('ch3sim_num/', alliters_sum[rep]))
  
  
  y_haj[,,,rep] <- clust_avg_y_sim[,,,1]
  wt_haj[,,,rep] <- hajek_wts[,,,1]
  num[,,,rep] <- sum_num[,,,1]
  
}
avg_y_haj = apply(y_haj, c(2,3), mean, na.rm = T)
lb_y_haj = apply(y_haj, c(2,3), quantile, probs = 0.025, na.rm = T)
ub_y_haj = apply(y_haj, c(2,3), quantile, probs = 0.975, na.rm = T)
apply(num, c(2,3), mean, na.rm = T)

avg_wt_haj = apply(wt_haj, c(2,3), mean, na.rm = T)

test = data.frame(g = rep(gamma_cts[2,], 2),
                  time = c(rep(1, n_gamma), rep(2, n_gamma)),
                  true_y = c(avg_truth[1,], avg_truth[2,]),
                  y_haj = c(avg_y_haj[1,], avg_y_haj[2,]),
                  lb_haj = c(lb_y_haj[1,], lb_y_haj[2,]),
                  ub_haj = c(ub_y_haj[1,], ub_y_haj[2,]))

ggplot(test, aes(x = g, y = y_haj)) + 
  geom_point() + 
  geom_line(aes(y = true_y)) +
  facet_wrap(~time) + 
  geom_ribbon(aes(ymin = lb_haj, ymax = ub_haj, xmin = g, xmax = g), alpha = 0.2, fill= 'blue')

#COVERAGE
truth_use = apply(truth, c(2,3,4), mean, na.rm = T) #avg across clusers

cov = c()
for(row in 1:nrow(test)){
  print(row)
  lb = test$lb_haj[row]
  ub = test$ub_haj[row]
  
  bw = c()
  for(i in 1:800){
    truth_vec = c(truth_use[1,,i], truth_use[2,,i])[row]
    bw = append(bw, between(truth_vec, lb, ub))
 }
cov = append(cov, mean(bw, na.rm = T))
}
cov
