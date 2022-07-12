#set working director
getwd()
setwd("~/project/cai")

#packages
library(devtools)
install_github("gpapadog/Interference", lib = getwd())

library(tidyverse)
library(Interference, lib.loc = getwd())
library(matrixStats)
library(abind)
library(parallel)

#scripts to load
source('GroupIPW_sd v2.R')
source('denominator_function_v2_sd.R')
source('helper_functs.R')
source('parallel_bootvar_function.R')
source('ypop_sd.R')
source('hajek_adj.R')
source('de_sd.R')
source('ie_sd.R')
source('CalcHeterogenousTrueIE.R')

#dataset to load
df = read.csv('cai_df_csv.csv')

#constants across simulations
estimand <- '1'
cov_cols = c(3)
ps_with_re = FALSE
numerator_with_re = FALSE
out_name = 'Yij'

phi_hat <- list(coefs = c(0, 0), re_var = 0)


beta_0 = .1
beta_1 = 3 #individual trt
beta_2 = 0 #covar
beta_3 = 0#.5 # proportion treated in cluster
p_trt = 0.5
giant_sim = TRUE


alpha <- c(0.5)

#get dims for simulations from OG data
if(giant_sim){
  test = df %>% group_by(neigh) %>% summarise(cluster_pop = n())
  test = rbind(test, test, test, test,test)
  test$neigh = c(1:nrow(test))
}else{
  test = df %>% group_by(neigh) %>% summarise(cluster_pop = n())
  test$neigh = c(1:nrow(test)) 
}

#test$cluster_pop = 15
#get possible gammas
gamma_list = rbind(0,seq(from = -.2, to = .2, length.out = 11))
gamma_numer = gamma_list

clust_size = '15'
get_yhats_boot <- function(clust_size){
  #simulate new data set
  sim_df = data.frame(neigh = c(), Aij = c(), Xij = c(), Tij = c())
  for (neigh in test$neigh){

    if(clust_size == 'obs'){village_size = test$cluster_pop[neigh]
    }else if(clust_size == '15'){village_size = 15#ifelse(neigh==1, 16, 15)#test$cluster_pop[neigh]
    }else if(clust_size == '50'){village_size = 50#ifelse(neigh==1, 51, 50)#test$cluster_pop[neigh]
    }

    #create Aij - the individual treatments
    Aij = rbinom(village_size, 1, p_trt) #this is trt vector for cluster
    Tij = sum(Aij) / village_size
    Xij = rnorm(village_size,0,1)
    XTij = Xij * Tij

    sim_df = rbind(sim_df, data.frame(neigh, Aij, Xij, Tij))
  }

  tst = beta_0 + beta_1*sim_df$Aij + beta_2*sim_df$Xij + beta_3*sim_df$Tij 
  sim_df$Yij = rnorm(nrow(sim_df), tst, 1)

  dta = sim_df
  cov_names = names(dta)[cov_cols]
  dta$Trt = dta$Aij
  
  #get the true heterogeneous interference results
  het_ie_truth = get_het_ie(dta, gamma_numer, cov_cols)
  

  trt_col <- which(names(dta) == 'Trt')
  out_col <- which(names(dta) == out_name)
  
  #get potential outcomes (not hajek)
  if(clust_size == 'obs'){
    yhat_group_sd2 <- GroupIPW_sd2(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                                   alpha = alpha, trt_col = trt_col, out_col = out_col,
                                   estimand = estimand, alpha_re_bound = 15, verbose = F,
                                   gamma_numer = gamma_list, loud_denom = F, fix_phi = T)
  }else{
    yhat_group_sd2 <- GroupIPW_sd2(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                                   alpha = alpha, trt_col = trt_col, out_col = out_col,
                                   estimand = estimand, alpha_re_bound = 20, verbose = F,
                                   gamma_numer = gamma_list, loud_denom = F, fix_phi = T, const_size = T)
  }

  z <- yhat_group_sd2$yhat_group
  po = apply(z, c(2,3), mean)

  haj = make_hajek(yhat_group_sd2, dta)

  #and bootstrap
  boots_est <- BootVar_sd(dta = dta,
                          B = 3,
                          alpha = alpha,
                          gamma_numer = gamma_list,
                          #ps = 'true',
                          cov_cols = cov_cols,
                          trt_col = trt_col,
                          out_col = out_col,
                          #ps_info_est = ps_info_est,
                          #phi_hat = phi_hat,
                          phi_hat_true = NULL,
                          verbose = TRUE,
                          return_everything =TRUE)
  
  

  ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = F, dta = dta)

  direct = DE_sd(ypop = ypop$ypop, ypop_var = ypop$ypop_var, boots = boots_est$hajboots)

  indirect0 = IE_sd(ygroup = yhat_group_sd2$yhat_group[,1,], boots = boots_est$hajboots, ps = 'true', hajofygroup = haj$haj[1,])[,,6]
  indirect1 = IE_sd(ygroup = yhat_group_sd2$yhat_group[,2,], boots = boots_est$hajboots, ps = 'true', hajofygroup = haj$haj[2,])[,,6]

  #horvitz thompson
  ht_ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = T, dta = dta)
  ht_direct = DE_sd(ypop = ht_ypop$ypop, ypop_var = ht_ypop$ypop_var, boots = boots_est$boots)

  #len = direct['boot_var_UB',] - direct['boot_var_LB',]
  coverage = (3 > direct['boot_var_LB',] &  3 < direct['boot_var_UB',])
  boot_se = sqrt(direct['boot_var',])

  coverage_ht = (3 > ht_direct['boot_var_LB',] &  3 < ht_direct['boot_var_UB',])
  boot_se_ht = sqrt(ht_direct['boot_var',])

  return(list(unadj = po,
              haj = haj$haj,
              wtlist = haj$wts,
              #interval = len,
              coverage=coverage,
              boot_se = boot_se,
              coverage_ht = coverage_ht,
              boot_se_ht = boot_se_ht,
              direct = direct,
              ht_direct = ht_direct,
              indirect0 = indirect0,
              indirect1 = indirect1))
}

#run above in parallel
n_sim <- 5
(n_cores <- detectCores() - 1)

output <- mclapply('15',mc.cores = n_cores, function(i) {
  res <- replicate(n_sim, get_yhats_boot(i))
  res#c(mean(res), sd(res))
}) #started at 12:48


#r = get_yhats_boot(clust_size = 'obs')
#r$wtlist

#get avg over many
n = 100
empty_unadj = list()#array(0, dim = c(2,11))
empty_haj = list()#array(0, dim = c(2,11))
empty_wts = list()

#empty_len = array(NA, dim = c(1,11,n))
empty_cov = array(NA, dim = c(1,11,n))
empty_bootse = array(NA, dim = c(1,11,n))

empty_cov_ht= array(NA, dim = c(1,11,n))
empty_bootse_ht = array(NA, dim = c(1,11,n))

#haj
empty_de = array(NA, dim=c(1,11,n))
empty_de_bootvar= array(0, dim=c(1,11,n))
empty_de_varform= array(0, dim=c(1,11,n))

empty_ie0 = array(NA, dim = c(1,11,n))
empty_ie0_lb = array(NA, dim = c(1,11,n))
empty_ie0_ub = array(NA, dim = c(1,11,n))

empty_ie1 = array(NA, dim = c(1,11,n))
empty_ie1_lb = array(NA, dim = c(1,11,n))
empty_ie1_ub = array(NA, dim = c(1,11,n))

#ht
empty_de_ht = array(NA, dim=c(1,11,n))
empty_de_bootvar_ht= array(0, dim=c(1,11,n))
empty_de_varform_ht= array(0, dim=c(1,11,n))

for(trial in 1:n){
  print(trial)
  t = get_yhats_boot(clust_size = '15')

  empty_unadj[[trial]] = t$unadj
  empty_haj[[trial]] = t$haj
  empty_wts[[trial]] = t$wtlist

  #hajek direct effect
  empty_de[,,trial] = t$direct[1,]
  empty_de_bootvar[,,trial] = t$direct[5,]
  empty_de_varform[,,trial] = t$direct[2,]

  #hajek indirect effect
  empty_ie0[,,trial] = t$indirect0[1,]
  empty_ie0_lb[,,trial] = t$indirect0['boot_var_LB',]
  empty_ie0_ub[,,trial] = t$indirect0['boot_var_UB',]

  empty_ie1[,,trial] = t$indirect1[1,]
  empty_ie1_lb[,,trial] = t$indirect1['boot_var_LB',]
  empty_ie1_ub[,,trial] = t$indirect1['boot_var_UB',]

  #ht direct effect
  empty_de_ht[,,trial] = t$ht_direct[1,]
  empty_de_bootvar_ht[,,trial] = t$ht_direct[5,]
  empty_de_varform_ht[,,trial] = t$ht_direct[2,]

  #empty_len <- rbind(empty_len, t$interval)
  empty_cov[,,trial] <-  t$coverage
  empty_bootse[,,trial] <- t$boot_se

  empty_cov_ht[,,trial] = t$coverage_ht
  empty_bootse_ht[,,trial] =t$boot_se_ht

  #if (n==300){save.image('feb18wkspc.Rdata')}
}
save.image('feb23wkspc.Rdata')
print(Sys.time()) #start at 3pm
load('feb18wkspc.Rdata')
hj_bias = apply(empty_de - 3, c(1,2), mean)
ht_bias = apply(empty_de_ht - 3, c(1,2), mean)

hj_cov = apply(empty_cov, c(2), mean)
ht_cov = apply(empty_cov_ht, c(2), mean)

hj_se = apply(empty_bootse, 2, mean)
ht_se = apply(empty_bootse_ht, 2, mean)

empty_haj = array(unlist(empty_haj), dim = c(2,10,n))
mc_var_haj = apply(empty_haj, c(1,2), var) #this is empirical se^2 estimated by monte carlo

empty_ie0 = apply(empty_ie0, c(1,2), mean)
empty_ie0_lb = apply(empty_ie0_lb, c(1,2), mean)
empty_ie0_ub = apply(empty_ie0_ub, c(1,2), mean)

empty_ie1 = apply(empty_ie1, c(1,2), mean)
empty_ie1_lb = apply(empty_ie1_lb, c(1,2), mean)
empty_ie1_ub = apply(empty_ie1_ub, c(1,2), mean)


ie0_df = data.frame(est = c(empty_ie0), lb = c(empty_ie0_lb), ub = c(empty_ie0_ub), gamma = gamma_list[2,], trt = '0')
ie1_df = data.frame(est = c(empty_ie1), lb = c(empty_ie1_lb), ub = c(empty_ie1_ub), gamma = gamma_list[2,], trt = '1')

ie_df = rbind(ie0_df, ie1_df)

ggplot(ie_df, aes(x = gamma, y = est, group = trt, col = trt)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  ggtitle('Indirect Effect') 
#ggsave('feb21/ie_188clusters_hajek_100sims_100boots.pdf')

#empty_unadj = array(unlist(empty_unadj), dim = c(2,10,n))
#mc_var_ht = apply(empty_unadj, c(1,2), var) #this is empirical se^2 estimated by monte carlo

#save.image('feb17wkspc.Rdata')
#load('feb18wkspc.Rdata')

hj_df <- data.frame(bias = c(hj_bias), coverage = hj_cov, boot_se = hj_se, gamma = gamma_list[2,], estimator = 'Hajek')
ht_df <- data.frame(bias = c(ht_bias), coverage = ht_cov, boot_se = ht_se, gamma = gamma_list[2,], estimator = 'Horvitz Thompson')

dd = rbind(hj_df, ht_df) %>% pivot_longer(cols = bias:boot_se)
ggplot(data = dd, 
       aes(x = gamma, y = value, color = estimator)) + 
  geom_point() + 
  facet_wrap(~name) + 
  theme(text = element_text(size = 14)) + 
  ggtitle('500 simulations, 100 bootstraps')
ggsave('feb21/bootstrapcovergae.pdf', width = 8, height = 5)


#compare mc se, boot se and formula se
#haj
mc_haj_var_de = apply(empty_de, c(1,2), var)
hj_var = apply(empty_de_varform, c(1,2), mean)
hj_bootvar = apply(empty_de_bootvar, c(1,2), mean)

#ht
mc_ht_var_de = apply(empty_de_ht, c(1,2), var)
ht_var = apply(empty_de_varform_ht, c(1,2), mean)
ht_bootvar = apply(empty_de_bootvar_ht, c(1,2), mean)

var_table_hj = data.frame(gamma = gamma_list[2,], 
                       mc_var = c(mc_haj_var_de),
                       formula_var = c(hj_var),
                       boot_var = c(hj_bootvar)) %>%
  pivot_longer(cols = mc_var:boot_var) %>%
  mutate(estimator = 'Hajek')

var_table_ht = data.frame(gamma = gamma_list[2,],
                          mc_var =c(mc_ht_var_de),
                          formula_var = c(ht_var),
                          boot_var = c(ht_bootvar)) %>%
  pivot_longer(cols = mc_var:boot_var) %>%
  mutate(estimator = 'HT')

var_table = rbind(var_table_hj, var_table_ht)

ggplot(data = var_table,
       aes(x = gamma, y = value, group = estimator, color = estimator)) + 
  facet_wrap(~name) + 
  geom_point() + 
  labs(title = 'Comparing Bootstrapped, Formulaic, \nand MC Variances of Hajek and HT DE',
       y = 'variance')
ggsave('feb21/variancecomp.pdf', width = 8, height = 5)

#look at just haj
ggplot(data = var_table %>% filter(estimator == 'Hajek'),
       aes(x = gamma, y = value, group = estimator, color = estimator)) + 
  facet_wrap(~name) + 
  geom_point() + 
  labs(title = 'Comparing Bootstrapped, Formulaic, \nand MC Variances of Hajek and HT DE',
       y = 'variance')
ggsave('feb21/variancecomp.pdf', width = 8, height = 5)

avg = array(0, dim = c(2,11))
avg2 = array(0, dim = c(2,11))
haj_de = array(0, dim = c(1,11))
ht_de = array(0, dim = c(1,11))
avg_wt = array(0, dim = c(2, 11))

#get mean over all and de dist
for(trial in 1:n){
  avg = avg + empty_haj[[trial]]
  avg2 = avg2 + empty_unadj[[trial]]
  
  haj_de = rbind(haj_de, empty_haj[[trial]][2,] - empty_haj[[trial]][1,])
  ht_de = rbind(ht_de, empty_unadj[[trial]][2,] - empty_unadj[[trial]][1,])
  
  avg_wt = avg_wt + empty_wts[[trial]]
}
avg = avg / n
avg2 = avg2/n
haj_de = haj_de[-1,]
ht_de = ht_de[-1,]

avg_wt = avg_wt / n

avg_df = data.frame(avg)
names(avg_df) = paste0('gamma=', gamma_list[2,])

avg2_df = data.frame(avg2)
names(avg2_df) = paste0('gamma=', gamma_list[2,])

avg_df = rbind(avg_df, avg2_df)

#plot the potential outcomes
#called avg (for haj adj)
#error bars
lb = array(0, dim = c(2,11))
ub = array(0, dim = c(2,11))
for (g in 1:11){
  for (z in 1:2){
    e = c()
    for (i in 1:n){
      e = append(e, empty_haj[[i]][z,g])
    }
    lb[z,g] = quantile(e, .025)
    ub[z,g] = quantile(e, .975)
  }
}

trt_df = data.frame(gamma = gamma_list[2,],
                    yhat = avg[2,], 
                    lower = lb[2,],
                    upper = ub[2,], trt = 'trt = 1')
untrt_df = data.frame(gamma = gamma_list[2,],
                      yhat = avg[1,], 
                      lower = lb[1,],
                      upper = ub[1,], trt = 'trt = 0')
p_df = rbind(trt_df, untrt_df)

ggplot(p_df, aes(x = gamma, y = yhat, col = trt)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~trt, scales = 'free_y') + 
  ggtitle('Avg Potential Outcomes - 3000 Simulations. \n47  clusters w cluster size 50. \nalpha = 0.5') 
ggsave('jan18_meeting/yhat_plot_3000_cluster50_error_homogint.pdf')
  
#plot the ie (should be flat bc no interference)
ie <- array(NA, dim = c(n, 2, 11, 11))
for (trial in 1:n){
  for(trt in 1:2){
    for (a1 in 1 : 11) {
      for (a2 in 1 : 11) {
        ie[trial, trt, a1, a2] <- empty_haj[[trial]][trt,a1] - empty_haj[[trial]][trt,a2]
      }
    }
  }
}

#untreated, 0 is reference gamma ie[,1,6,]
ie0 = data.frame(t(apply(ie[,1,6,], 2, quantile, c(.025, .5, .975))))
names(ie0) = c('lb', 'est', 'ub')
ie0$gamma = gamma_list[2,]; ie0$trt = '0'
#treated, 0 is reference gamma ie[,2,6,]
ie1 = data.frame(t(apply(ie[,2,6,], 2, quantile, c(.025, .5, .975))))
names(ie1) = names(ie0)[1:3]
ie1$gamma = gamma_list[2,]; ie1$trt = '1'
ie = rbind(ie0, ie1)

ggplot(ie, aes(x = gamma, y = est, col = trt)) + 
  geom_point() + 
  #geom_errorbar(aes(ymin = lb, ymax = ub)) +
  facet_wrap(~trt, scales = 'free_y') + 
  ggtitle('IE - 3000 Simulations. \n47  clusters w cluster size 50. \nalpha = 0.5') 
ggsave('jan18_meeting/ie_plot_3000_cluster50.pdf')



#plot the de
de = apply(haj_de, 2, mean)
de_lb = colQuantiles(haj_de, probs = 0.025)
de_ub = colQuantiles(haj_de, probs = 0.975)
de = data.frame(de, de_lb, de_ub) %>%
  mutate(gamma = gamma_list[2,])

ggplot(data = de, aes(x = gamma, y = de)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = de_lb, ymax = de_ub)) + 
  ggtitle('Direct Effect. Simulate n = 3000 \n47 simulated clusters w cluster size 50. \nalpha = 0.5') + 
  ylab('DE')
ggsave('jan5_meeting/de_plot_3000_cluster50.pdf')

#plot Horvitz Thompson DE
#plot the de
de_ht = apply(ht_de, 2, mean)
de_lb_ht = colQuantiles(ht_de, probs = 0.025)
de_ub_ht = colQuantiles(ht_de, probs = 0.975)
de_ht = data.frame(de_ht, de_lb_ht, de_ub_ht) %>%
  mutate(gamma = gamma_list[2,])

ggplot(data = de_ht, aes(x = gamma, y = de_ht)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = de_lb_ht, ymax = de_ub_ht)) + 
  ggtitle('HT Direct Effect. Simulate n = 3000 \n47 simulated clusters w cluster size 50. \nalpha = 0.5') + 
  ylab('DE')
ggsave('jan5_meeting/ht_de_plot_3000_cluster50.pdf')




#a boxplot for each gamma / trt combo where i have a point for each sim







#plot the sum of weights thing GA wants
size = 50*47#sum(test$cluster_pop)
plot_df = data.frame()
for (g in 1:11){
  for (z in 1:2){
    wts = avg_wt[z,g]
    wtdf = data.frame(neigh = size, weightsum = wts) %>% 
      mutate(lab = 'Sum Weights', gamma = gamma_list[2,g], z = z-1)
    plot_df = rbind(plot_df, wtdf)
  }
}

#error bars
lb = array(0, dim = c(2,11))
ub = array(0, dim = c(2,11))
for (g in 1:11){
  for (z in 1:2){
    e = c()
    for (i in 1:n){
      e = append(e, empty_wts[[i]][z,g])
    }
    lb[z,g] = quantile(e, .025)
    ub[z,g] = quantile(e, .975)
  }
}


lb_df = rbind(data.frame(gamma = gamma_list[2,], 
                         lb = lb[1,]/size, 
                         z = 0),
      data.frame(gamma = gamma_list[2,], 
                 lb = lb[2,]/size, 
                 z = 1))

ub_df = rbind(data.frame(gamma = gamma_list[2,], ub = (ub[1,]/size), z = 0),
              data.frame(gamma = gamma_list[2,], ub = ub[2,]/size, z = 1))

plot_df = plot_df %>% mutate(wt_pct = weightsum / size) %>%
  merge(lb_df, by = c('gamma', 'z')) %>%
  merge(ub_df, by = c('gamma', 'z'))

ggplot(data = plot_df,
       aes(x = gamma, y = wt_pct)) + 
  geom_point(alpha = 0.5) +
  #geom_errorbar(aes(ymin = lb, ymax = ub)) + 
  facet_wrap(~z) + 
  geom_hline(yintercept = 1) + 
  ggtitle('Sum of weights as percent of population size (n = 3000) \ncluster size 50') + 
  ylab('sum of weights / pop size')
ggsave('jan5_meeting/wts_by_popsize_3000_clustersize50.pdf')


save.image('Jan5wkspc.Rdata')


#does the shape of de look the same for different trials?
plot(haj_de[1,], type = 'l', ylab = 'DE', xlab = 'Gamma', main = "Does the shape of the DE curve mean anything?", ylim = c(2.85, 3.1))
for (i in 1:20){
  points(haj_de[i,], type = 'l', col = i)
}
ggsave('jan5_meeting/DEshape_3000_15clustersize.pdf')










#plot boxplot of estimators by gamma, trt for the 500 trials
plot_df_box = data.frame(matrix(ncol = 3, nrow = 0))
names(plot_df_box) = c('gamma', 'trt', 'yhat')
for (trial in 1:n){
  for (g in 1:11){
    for (z in 1:2){
      target = empty_haj[[trial]][z,g]
      
      dfrow = data.frame(gamma =g, trt = z, yhat = target)
      plot_df_box = rbind(plot_df_box, dfrow)

    }
  }
}

ggplot(data = plot_df_box, aes(y = yhat)) + 
  geom_boxplot() + 
  facet_grid(trt~gamma, scales = 'free_y')
ggsave('boxplot of yhats_5000.pdf', width = 8, height = 6)
