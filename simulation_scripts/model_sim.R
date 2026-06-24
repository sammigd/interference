print(getwd())

bivar = F
alt_clus_size = F

package_path <- "/nfs/roberts/project/pi_lf474/sgd37/pckgs"  # adjust as needed

#packages
library(usethis)
library(rlang)
library(devtools)
#install_github("gpapadog/Interference", lib = package_path)
library(tidyverse)
library(Interference, lib.loc = package_path)
library(matrixStats)
library(abind)
library(parallel)
options(dplyr.summarise.inform = FALSE)


#set working director
getwd()
#setwd("~/project/cai")
setwd('/nfs/roberts/project/pi_lf474/sgd37')

#scripts to load
source('interference/analysis_scripts/GroupIPW_sd v2.R')
source('interference/analysis_scripts/denominator_sd.R')
source('interference/analysis_scripts/helper_functs.R')
source('interference/analysis_scripts/parallel_bootvar_function test.R')
source('interference/analysis_scripts/ypop_sd.R')
source('interference/analysis_scripts/hajek_adj.R')
source('interference/analysis_scripts/de_sd.R')
source('interference/analysis_scripts/ie_sd.R')
source('interference/analysis_scripts/CalcHeterogeneousTrueIE.R')
source('interference/analysis_scripts/oe_sd.R')
#source('scripts/diffusion_truth.R')

print('packages and scripts loaded')

######################################################
# PARAMETER SET UP ###################################
######################################################

#get args from slurm
args <- commandArgs(TRUE)
#args = c(1, 0, 0, 2, 0, 1)
#args are 1: taskid
#         2: diffusion
#         3: b3
#         4: b4
#         5: concordance
#         6: b5
print(args)
idx = as.numeric(args[1])
diffusion = as.logical(as.numeric(args[2]))
if(diffusion){bivar = F}
  

if(!diffusion){
  #concor = 0; beta_3 = 0; beta_4 = 0; beta_5 = 0
  beta_3 = as.numeric(args[3]) #0,1
  beta_4 = as.numeric(args[4]) #0,.5, 1
  concor = as.numeric(args[5]) #0, .65, .8
  beta_5 = as.numeric(args[6]) #for bivariate X2 spillover
  hypothetical_alpha = as.numeric(args[7])
  ngamma = as.numeric(args[8])
  print(paste('ngamma = ', ngamma))
  
  clust_size = '15'
  n_clus = 200
  epsilon = 1

}else{ #diffusion == T
  concor = as.numeric(args[3])
  diffusion_p = as.numeric(args[4])
  clust_size = 5
  n_clus = 200
}



print('setup done')
#dataset to load
#df = read.csv('cai_df_csv.csv')


###################################
# SETUP CONSTANTS #################
###################################


#constants across simulations
estimand <- '1'
if(!diffusion){
  cov_cols = c(3,4,5)
  beta_0 = .1
  beta_1 = 3 #individual trt
  beta_2 = 0 #covar
  observed_alpha = 0.5
  #hypothetical_alpha = 0.35 #0.65
}else{ #if diffusion
  cov_cols = c(3,4)
  observed_alpha = 0.25
  hypothetical_alpha = 0.25
  x2_prev = 0.2
}

out_name = 'Yij'

test = data.frame(neigh = c(1:n_clus),
                  cluster_pop = rep(15, n_clus))
beta_6 = 0
if(alt_clus_size == T){
  test = data.frame(neigh = c(1:n_clus),
                    cluster_pop = c(rep(15, n_clus/2), c(rep(10, n_clus/2))))
  beta_6 = 2
}

#get possible gammas
if(!diffusion){
  if(!bivar){ #univar
    gl = seq(from = -1.3, to = 1.3, length.out = ngamma)
    ngl = rep(0, ngamma)
    gamma_list = rbind(rep(0, ngamma*3),
                      c(gl, ngl, ngl),
                      c(ngl, gl, ngl),
                      c(ngl, ngl, gl))
    gamma_list = cbind(gamma_list, c(0,0,0,0))
    
  }
  if(bivar){ #bivar
    gl = c(seq(from = -1.3, to = 1.3, length.out = 7), 0)
    ggrid = expand.grid(gl, gl)
    gamma_list = t(cbind(0, ggrid, 0))
  }
  
  
  
}
if(diffusion){
  gl = seq(from = -.2, to = .2, length.out = 33)
  ngl = rep(0, 33)
  gamma_list = rbind(rep(0, 33),
                     c(ngl, gl),
                     c(gl, ngl))
  gamma_list = cbind(gamma_list, c(0,0,0))
}

gamma_numer = gamma_list
ngam = ncol(gamma_numer)
n_boot = 2#150

#plot(gamma_numer[2,], gamma_numer[3,])


##########################################################
# SETTING SEED ###########################################
##########################################################
#idx is 1:600 for each set of parms
if(diffusion){
  print('is diffusion')
  if(concor == 0 & diffusion_p == 0){c_seed = 1} 
  if(concor == 0 & diffusion_p == 0.2){c_seed = 1} 
  if(concor == 0 & diffusion_p == 0.5){c_seed = 3} 
  if(concor == 0 & diffusion_p == 0.8){c_seed = 4} 
  if(concor == 0.65 & diffusion_p == 0){c_seed = 5}
  if(concor == 0.65 & diffusion_p == 0.2){c_seed = 6}
  if(concor == 0.65 & diffusion_p == 0.5){c_seed = 7}
  if(concor == 0.65 & diffusion_p == 0.8){c_seed = 8}
}
if(!diffusion){
  print('not diffusion')
  if(concor == 0 & beta_3 == 0 & beta_4 == 0 & beta_5 == 0){c_seed = 9} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 0 & beta_5 == 0){c_seed = 10} 
  if(concor == 0 & beta_3 == 0 & beta_4 == 1 & beta_5 == 0){c_seed = 11} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 1 & beta_5 == 0){c_seed = 12} 
  if(concor == 0 & beta_3 == 0 & beta_4 == 2 & beta_5 == 0){c_seed = 13} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 2 & beta_5 == 0){c_seed = 14} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 0 & beta_5 == 0){c_seed = 15} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 0 & beta_5 == 0){c_seed = 16} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 1 & beta_5 == 0){c_seed = 17} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 1 & beta_5 == 0){c_seed = 18} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 2 & beta_5 == 0){c_seed = 19} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 2 & beta_5 == 0){c_seed = 20} 
  if(concor == 0 & beta_3 == 0 & beta_4 == 0 & beta_5 == 1){c_seed = 21} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 0 & beta_5 == 1){c_seed = 22} 
  if(concor == 0 & beta_3 == 0 & beta_4 == 1 & beta_5 == 1){c_seed = 23} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 1 & beta_5 == 1){c_seed = 24} 
  if(concor == 0 & beta_3 == 0 & beta_4 == 2 & beta_5 == 1){c_seed = 25} 
  if(concor == 0 & beta_3 == 1 & beta_4 == 2 & beta_5 == 1){c_seed = 26} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 0 & beta_5 == 1){c_seed = 27} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 0 & beta_5 == 1){c_seed = 28} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 1 & beta_5 == 1){c_seed = 29} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 1 & beta_5 == 1){c_seed = 30} 
  if(concor == 0.65 & beta_3 == 0 & beta_4 == 2 & beta_5 == 1){c_seed = 31} 
  if(concor == 0.65 & beta_3 == 1 & beta_4 == 2 & beta_5 == 1){c_seed = 32} 
}

print(c_seed)

get_yhats_boot <- function(clust_size, 
                           use.boot = T, 
                           n_boot, 
                           epsilon = NA, 
                           diffusion, 
                           nn, 
                           diffusion_p=NULL, 
                           seed = NULL){
  #simulate new data set
  set.seed(seed)
  if(!diffusion){
    sim_df = data.frame(neigh = c(), 
                        Aij = c(), 
                        X1ij = c(), X2ij = c(), X3ij = c(), 
                        Tij = c(), Tprimeij =c(), Tprimeij2 =c())
    for (neigh in test$neigh){
      #neigh = 1
      if(clust_size == 'obs'){village_size = test$cluster_pop[neigh]
      }else if(clust_size == '15'){village_size = 15#ifelse(neigh==1, 16, 15)#test$cluster_pop[neigh]
      }else if(clust_size == '5'){village_size = 5#ifelse(neigh==1, 51, 50)#test$cluster_pop[neigh]
      }
      
      if(alt_clus_size == T){
        village_size = test$cluster_pop[neigh]
      }
      
      #create Aij - the individual treatments
      Aij = rbinom(village_size, 1, observed_alpha) #this is trt vector for cluster
      if(sum(Aij) == 0 | sum(Aij) == length(Aij)){Aij = rbinom(village_size, 1, observed_alpha)}
      X1ij = rbinom(village_size,1, 0.5)
      X3ij = rbinom(village_size,1, 0.5)
      
      
      X2ij = rbinom(village_size,1, 0.5)
      if(concor == 0.65){ 
        X2ij[X1ij == 1 & X2ij == 0] <- rbinom(length(X2ij[X1ij == 1 & X2ij == 0]), 1, .255)
        X2ij[X1ij == 0 & X2ij == 1] <- rbinom(length(X2ij[X1ij == 0 & X2ij == 1]), 1, .745)
      }
      if(concor == 0.8){
        X2ij[X1ij == 1 & X2ij == 0] <- rbinom(length(X2ij[X1ij == 1 & X2ij == 0]), 1, .55)
        X2ij[X1ij == 0 & X2ij == 1] <- rbinom(length(X2ij[X1ij == 0 & X2ij == 1]), 1, .45)
      }
      
      Tij = (sum(Aij) - Aij) / (village_size-1)
      Tprimeij = (sum(Aij*X1ij) - (Aij*X1ij)) / (sum(Aij)-Aij)
      Tprimeij = ifelse(sum(Aij) - Aij == 0, 0, Tprimeij)
      
      Tprimeij2 = (sum(Aij*X2ij) - (Aij*X2ij)) / (sum(Aij)-Aij)
      Tprimeij2 = ifelse(sum(Aij) - Aij == 0, 0, Tprimeij2)
      
      sim_df = rbind(sim_df, 
                     data.frame(neigh = rep(neigh, village_size), 
                                Aij, 
                                X1ij, X2ij, X3ij, 
                                Tij, Tprimeij, Tprimeij2))
    }

    tst = beta_0 + beta_1*sim_df$Aij + beta_2*sim_df$X1ij + beta_3*sim_df$Tij + beta_4*sim_df$Tprimeij + beta_5*sim_df$Tprimeij2 + beta_6*ifelse(village_size == 15, 1, 0)*sim_df$Aij
    sim_df$Yij = rnorm(nrow(sim_df), tst, epsilon)
    
    dta <- sim_df
  }
  if(diffusion){
    neigh = sort(rep(1:n_clus, clust_size))
    X1ij = rep(c(1,0,0,0,0), n_clus)
    Aij = rbinom(n_clus*clust_size, 1, observed_alpha)
    Yij = rep(NA, n_clus*clust_size)
    
    sim = data.frame(neigh, X1ij, Aij, Yij)
    
    #gen outcome from analytical probs
    count_untrt_center = length(sim$Yij[sim$Aij == 0 & sim$X1ij == 1])
    count_untrt_alter = length(sim$Yij[sim$Aij == 0 & sim$X1ij == 0])
    
    
    sim$Yij[sim$Aij == 1] <- 1
    sim$Yij[sim$Aij == 0 & sim$X1ij == 1] <- rbinom(count_untrt_center, 1, 0.413)
    sim$Yij[sim$Aij == 0 & sim$X1ij == 0] <- rbinom(count_untrt_alter, 1, 0.125)
    
    sim$X2ij = 0
    
    #generate outcome
    sim$is_center_trted <- rep(sim %>% group_by(neigh) %>% slice(1) %>% pull(Aij), each = 5)
    sim$Tij <- rep(sim %>% group_by(neigh) %>% summarise(Tij = sum(Aij)) %>% pull(Tij), each = 5)
    
    for (i in 1:nrow(sim)){
      if(sim$Aij[i] == 1){sim$Yij[i] = 1}
      
      #untreated alter
      if(sim$Aij[i] == 0 & sim$X1ij[i] == 0 & sim$is_center_trted[i] == 0){sim$Yij[i] = 0}
      if(sim$Aij[i] == 0 & sim$X1ij[i] == 0 & sim$is_center_trted[i] == 1){sim$Yij[i] = rbinom(1,1,diffusion_p)}
      
      #untreated center
      if(sim$Aij[i] == 0 & sim$X1ij[i] == 1 & sim$Tij[i] == 0){sim$Yij[i] = 0}
      if(sim$Aij[i] == 0 & sim$X1ij[i] == 1 & sim$Tij[i] > 0){sim$Yij[i] = max(rbinom(sim$Tij[i],1,diffusion_p))}
    }
      
    #generate x2
    sim$X2ij = rbinom(clust_size * n_clus, 1, x2_prev) #kappa = 0
    
    if(concor == 0.65){
      kappa = 0.25
      sim$X2ij = rbinom(n_clus*clust_size, 1, 
                        rep(c(x2_prev + kappa, 
                              rep(x2_prev - kappa / (clust_size-1), clust_size-1)), 
                            n_clus)) 
    }
    
    dta <- sim %>% select(neigh, Aij, X1ij, X2ij, Yij, Tij)
      }
  cov_names  = names(dta)[cov_cols]
  dta$Trt = dta$Aij
  
  #DONE WITH DATA GENERATION - NOW RUN CODE
  
  #true_ygroup = dta %>% group_by(neigh) %>% summarise(Ybar = mean(Yij))
  #true_ygroup = matrix(true_ygroup$Ybar, ncol = 1)
  
  #get the true heterogeneous interference results
  print('start het truth')
  if(!diffusion){
    truth = get_het_ie(dta, gamma_numer, cov_cols,
                       beta_0 = beta_0, beta_1 = beta_1, beta_2 = beta_2, beta_3 = beta_3, beta_4 = beta_4, beta_5 = beta_5,
                       alpha = hypothetical_alpha, diffusion = diffusion, nn = nn)}
  if(diffusion){
    truth = get_het_ie(dta, gamma_numer, cov_cols, 
                       alpha = hypothetical_alpha, diffusion = diffusion, diffusion_p = diffusion_p, nn = nn)}#NULL}

  trt_col <- which(names(dta) == 'Trt')
  out_col <- which(names(dta) == out_name)
  
  #get potential outcomes (not hajek)
  print('start yhat')
  if(clust_size == 'obs'){
    print('yhats')
    yhat_group_sd2 <- GroupIPW_sd2(dta = dta, cov_cols = cov_cols,
                                   alpha = hypothetical_alpha, trt_col = trt_col, out_col = out_col,
                                   alpha_re_bound = 15, verbose = F,
                                   gamma_numer = gamma_list, loud_denom = F, fix_phi = T)
  }else{
    yhat_group_sd2 <- GroupIPW_sd2(dta = dta, cov_cols = cov_cols,
                                   alpha = hypothetical_alpha, trt_col = trt_col, out_col = out_col,
                                   alpha_re_bound = 20, verbose = F,
                                   gamma_numer = gamma_list, loud_denom = F, fix_phi = T, const_size = T)
  }

  z <- yhat_group_sd2$yhat_group
  po = apply(z, c(2,3), mean)
  
  #print('haj')
  haj = make_hajek(yhat_group_sd2, dta)

  #and bootstrap
  if(use.boot == T){
    print('start boot') #n_boot = 3
    boots_est <- BootVar_sd(dta = dta, B = n_boot, alpha = hypothetical_alpha, gamma_numer = gamma_list,
                            cov_cols = cov_cols, trt_col = trt_col, out_col = out_col,
                            phi_hat_true = NULL, verbose = TRUE, return_everything = TRUE)
    print('boot done')
    
  }else{boots_est = NULL}
  
  print('ypop')
  ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = F, dta = dta)
  #ypop bootvar: bootvariance for each ypop by gamma
  ypop_bootvar = apply(boots_est$boots, c(1,2), var)


  print('de')
  direct = DE_sd(ypop = ypop$ypop, ypop_var = ypop$ypop_var, boots = boots_est$hajboots, ygroup = yhat_group_sd2$yhat_group)
  direct_cov = ypop$ypop_var_de
  
  print('ie')
  indirect0 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,1,], 
                    hajofygroup = haj$haj[1,], boots = boots_est$hajboots, 
                    treatment = 1)$ie[,,ngam]
  indirect1 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,2,], 
                    hajofygroup = haj$haj[2,], boots = boots_est$hajboots, treatment =2)$ie[,,ngam]
  
  print('oe')
  oe = OE_sd(ypop = ypop, ygroup = yhat_group_sd2$oe_yhat_group,  
             hajofygroup = haj$oe_haj, boots = boots_est$oehajboots)
  oe_cov = oe$oe_var
  oe = oe$oe[,,ngam]

  #horvitz thompson
  print('htypop')
  ht_ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = T, dta = dta)
  
  print('htdirect')
  ht_direct = DE_sd(ypop = ht_ypop$ypop, ypop_var = ht_ypop$ypop_var, boots = boots_est$boots, ygroup = yhat_group_sd2$yhat_group)
  
  print('htindirect')
  ht_indirect0 =  IE_sd(ypop = ht_ypop$ypop[1,], ygroup = yhat_group_sd2$yhat_group[,1,], 
                        hajofygroup = haj$haj[1,], boots = boots_est$boots, horvitzthompson = T, treatment = 1)$ie[,,ngam]
  ht_indirect1 =  IE_sd(ypop = ht_ypop$ypop[2,], ygroup = yhat_group_sd2$yhat_group[,2,], 
                        hajofygroup = haj$haj[2,], boots = boots_est$boots, horvitzthompson = T, treatment = 2)$ie[,,ngam]

  ht_oe = OE_sd(ypop = ht_ypop, ygroup = yhat_group_sd2$oe_yhat_group, 
                hajofygroup = haj$oe_haj, boots = boots_est$oehtboots, horvitzthompson = T)
  ht_oe_cov = ht_oe$oe_var
  ht_oe = ht_oe$oe[,,ngam]
  
  ht_oe_mcvar = apply(boots_est$oe_est_ht, 1, var)
  haj_oe_mcvar = apply(boots_est$oe_est_ht,1, var)
  
  if(use.boot==T){
    
    return(list(unadj = po,
                haj = haj$haj,#this is y0 and y1
                oe_haj = haj$oe_haj, #this is y_avg
                wtlist = haj$wts,
                direct = direct,
                direct_cov = direct_cov,
                indirect_cov = ypop$ypop_var_ie,
                overall_cov = ypop$ypop_var_oe,
                ht_direct = ht_direct,
                indirect0 = indirect0,
                indirect1 = indirect1,
                ht_indirect0=ht_indirect0,
                ht_indirect1=ht_indirect1,
                ht_oe=ht_oe, ht_oe_cov = ht_oe_cov,
                truth = truth,
                oe = oe, oe_cov = oe_cov,
                ht_oe_mcvar = ht_oe_mcvar,
                haj_oe_mcvar = haj_oe_mcvar,
                ht_yhat_var = ht_ypop$ypop_var,
                ypop_bootvar = ypop_bootvar))
  }
  

  return(list(unadj = po,
              haj = haj$haj,
              wtlist = haj$wts,
              direct = direct,
              ht_direct = ht_direct,
              indirect0 = indirect0,
              indirect1 = indirect1,
              het_ie_truth = het_ie_truth))
}


if(diffusion){
  test = get_yhats_boot(clust_size, 
                        use.boot = T, 
                        n_boot = n_boot, 
                        epsilon = epsilon, 
                        diffusion = diffusion, 
                        nn = 100, 
                        diffusion_p = diffusion_p, 
                        seed = as.numeric(paste0(c_seed, 0, idx)))
  save.image(paste0('ms_diffusion_2024may11', '/scenario', idx, '_', concor,'_',diffusion_p, '.RSave'))
}else{
  test = get_yhats_boot(clust_size, 
                        use.boot = T, 
                        n_boot = n_boot, 
                        epsilon = epsilon, 
                        diffusion = diffusion, 
                        nn = 2,#100,
                        seed = as.numeric(paste0(c_seed, 0, idx)))
  if(bivar){save.image(paste0('ms_bivar_2024jun01', '/scenario', idx, '_',beta_3, '_', beta_4,'_', concor,'_', beta_5, '.RSave'))}
  if(!bivar){
    if(hypothetical_alpha == 0.4){
      save.image(paste0('ms_univar_2025jul9_smallalpha','/scenario', idx, '_',beta_3, '_', beta_4,'_', concor,'_', beta_5, '.RSave'))
    }
    if(hypothetical_alpha == 0.6){
      save.image(paste0('ms_univar_2025jul9_largealpha','/scenario', idx, '_',beta_3, '_', beta_4,'_', concor,'_', beta_5, '.RSave'))
    }
    if(hypothetical_alpha == 0.5){
      print('Made it to save step!')
      print('Working Directory:')
      print(getwd())
      if(alt_clus_size == F){save.image(paste0('interference/results/ms_univar_2026may16_stdalpha','/scenario', idx, '_',beta_3, '_', beta_4,'_', concor,'_', beta_5, '_', ngam, '.RSave'))}
      if(alt_clus_size == T){save.image(paste0('ms_univar_2025aug17_altclus','/scenario', idx, '_',beta_3, '_', beta_4,'_', concor,'_', beta_5, '.RSave'))}
      print('Saved!')
    }
  }
}



