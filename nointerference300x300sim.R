#set working director
getwd()
setwd("~/project/cai")

#packages
#install_github("gpapadog/Interference")
library(devtools)
library(tidyverse)
library(Interference)
library(matrixStats)

#scripts to load
source('GroupIPW_sd v2.R')
source('denominator_function_v2_sd.R')
source('helper_functs.R')
source('BootVar_sd.R')
source('ypop_sd.R')
source('hajek_adj.R')
source('de_sd.R')
source('ie_sd.R')
source('CalcHeterogenousTrueIE.R')

#dataset to load
df = read.csv('cai_df_csv.csv')

#constants across simulations
estimand <- '1'
cov_cols = c(3,4,5)
ps_with_re = FALSE
numerator_with_re = FALSE
out_name = 'Yij'

phi_hat <- list(coefs = c(0, 0), re_var = 0)


beta_0 = .1
beta_1 = 3 #individual trt
beta_2 = 0 #covar
beta_3 = 0 #.5 # proportion treated in cluster
beta_4 = 0#1 # interaction bw trt and group trt #heterogeneous interference
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
n_boot = 300

get_yhats_boot <- function(clust_size){
  #simulate new data set
  sim_df = data.frame(neigh = c(), Aij = c(), Xij = c(), Tij = c(), XTij = c())
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
    
    sim_df = rbind(sim_df, data.frame(neigh, Aij, Xij, Tij, XTij))
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
                          B = n_boot,
                          alpha = alpha,
                          gamma_numer = gamma_list,
                          #ps = 'true',
                          cov_cols = cov_cols,
                          trt_col = trt_col,
                          out_col = out_col,
                          #ps_info_est = ps_info_est,
                          phi_hat = phi_hat,
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
n_sim <- 300
(n_cores <- detectCores() - 1)

output <- pbmclapply('15',mc.cores = n_cores, function(i) {
  res <- replicate(n_sim, get_yhats_boot(i))
  res#c(mean(res), sd(res))
}) #started at 12:48

save.image('nointerference.RData')
