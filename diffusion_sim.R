#get args from slurm
args <- commandArgs(TRUE)
print(args)

idx = as.numeric(args[1])
concor = as.numeric(args[2]) #0, .65, .8
#concor = 0

#set working director
getwd()
setwd("~/project/cai")


#packages
library(devtools)
#install_github("gpapadog/Interference", lib = getwd())
library(tidyverse)
library(Interference, lib.loc = getwd())
library(matrixStats)
library(abind)
library(parallel)

#scripts to load
source('scripts/GroupIPW_sd v2.R')
source('scripts/denominator_sd.R')
source('scripts/helper_functs.R')
source('scripts/parallel_bootvar_function test.R')
source('scripts/ypop_sd.R')
source('scripts/hajek_adj.R')
source('scripts/de_sd.R')
source('scripts/ie_sd.R')
source('scripts/CalcHeterogenousTrueIE.R')
source('scripts/oe_sd.R')


#constants across simulations
estimand <- '1'
cov_cols = c(3,4)#,4,5)
ps_with_re = FALSE
numerator_with_re = FALSE
out_name = 'Yij'

phi_hat <- list(coefs = c(0, 0), re_var = 0)

p_trt = 0.25
alpha <- c(0.25)
diffusion_p = 0.5
x2_prev = 0.2

clust_size = 5
n_boot = 3#150
n_clus = 200

#gammas
gl = seq(from = -.2, to = .2, length.out = 33)
ngl = rep(0, 33)
gamma_list = rbind(rep(0, 33),
                   c(ngl, gl),
                   c(gl, ngl))
gamma_list = cbind(gamma_list, c(0,0,0))
gamma_numer = gamma_list

get_yhats_boot <- function(clust_size, use.boot = T){
  #generate new data set
  neigh = sort(rep(1:n_clus, clust_size))
  X1ij = rep(c(1,0,0,0,0), 200)
  Aij = rbinom(1000, 1, p_trt)
  Yij = rep(NA, 1000)
  
  sim = data.frame(neigh, X1ij, Aij, Yij)
  
  #generate outcome
  sim$is_center_trted <- rep(sim |> group_by(neigh) |> slice(1) |> pull(Aij), each = 5)
  sim$Tij <- rep(sim |> group_by(neigh) |> summarise(Tij = sum(Aij)) |> pull(Tij), each = 5)
  
  sim = sim %>%
    group_by(neigh) |>
    mutate(Yij = case_when(Aij == 1 ~ 1, #treated node
                           #untrt alter where center IS NOT treated
                           Aij == 0 & X1ij == 0 & is_center_trted == 0 ~ 0,
                           #untreated center - #max of a coin flip for each treated alter in the cluster
                           Aij == 0 & X1ij == 1 ~ as.numeric(max(rbinom(Tij, 1, diffusion_p))), 
                           #untrt alter where center IS treated
                           Aij == 0 & X1ij == 0 & is_center_trted == 1 ~ as.numeric(rbinom(1,1,diffusion_p))))
  
  sim = sim %>% group_by(neigh) %>%
    mutate(Tprimeij = sum(Aij[X1ij == 1]) / sum(Aij))
  
  #generate x2
  sim$X2ij = rbinom(1000, 1, x2_prev) #kappa = 0
  
  if(concor == 0.65){
    kappa = 0.25
    sim$X2ij = rbinom(n_clus*clust_size, 1, 
                  rep(c(x2_prev + kappa, 
                        rep(x2_prev - kappa / (clust_size-1), clust_size-1)), 
                      n_clus)) 
  }
  
  #for validating cohens kapps
  #kaps = c()
  #for (cluster in 1:200){
  #  kaps = append(kaps, calc_cohen_k(sim$X1ij[sim$neigh == cluster], sim$X2ij[sim$neigh == cluster]))
  #}
  #mean(kaps)
  
  dta <- sim %>% select(neigh, Aij, X1ij, X2ij, Yij, Tij, Tprimeij)
  cov_names  = names(dta)[cov_cols]
  dta$Trt = dta$Aij
  
  #true_ygroup = dta %>% group_by(neigh) %>% summarise(Ybar = mean(Yij))
  #true_ygroup = matrix(true_ygroup$Ybar, ncol = 1)
  
  #get the true heterogeneous interference results
  print('start het truth')
  truth = get_true_diffusion(gamma_numer = gamma_numer, 
                                    x2_prev = 0.2, 
                                    kappa = 0.25)
  #het_ie_truth = NULL
  
  trt_col <- which(names(dta) == 'Trt')
  out_col <- which(names(dta) == out_name)
  
  #get potential outcomes (not hajek)
  print('start yhat')
  if(clust_size == 'obs'){
    print('yhats')
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
  
  #print('haj')
  haj = make_hajek(yhat_group_sd2, dta)
  
  #and bootstrap
  if(use.boot == T){
    print('start boot') #n_boot = 3
    boots_est <- BootVar_sd(dta = dta, B = n_boot, alpha = alpha, gamma_numer = gamma_list,
                            cov_cols = cov_cols, trt_col = trt_col, out_col = out_col,
                            phi_hat_true = NULL, verbose = TRUE, return_everything = TRUE)
    print('boot done')
    
  }else{boots_est = NULL}
  
  print('ypop')
  ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = F, dta = dta)
  
  print('de')
  direct = DE_sd(ypop = ypop$ypop, ypop_var = ypop$ypop_var, boots = boots_est$hajboots)
  
  print('ie')
  indirect0 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,1,], 
                    hajofygroup = haj$haj[1,], boots = boots_est$hajboots, 
                    treatment = 1)[,,ngam]
  indirect1 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,2,], 
                    hajofygroup = haj$haj[2,], boots = boots_est$hajboots, treatment =2)[,,ngam]
  
  print('oe')
  oe = OE_sd(ypop = ypop, ygroup = yhat_group_sd2$oe_yhat_group,  
             hajofygroup = haj$oe_haj, boots = boots_est$oehajboots)
  oe_cov = oe$oe_var
  oe = oe$oe[,,ngam]
  
  #horvitz thompson
  print('htypop')
  ht_ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = T, dta = dta)
  
  print('htdirect')
  ht_direct = DE_sd(ypop = ht_ypop$ypop, ypop_var = ht_ypop$ypop_var, boots = boots_est$boots)
  
  print('htindirect')
  ht_indirect0 =  IE_sd(ypop = ht_ypop$ypop[1,], ygroup = yhat_group_sd2$yhat_group[,1,], 
                        hajofygroup = haj$haj[1,], boots = boots_est$boots, horvitzthompson = T, treatment = 1)[,,ngam]
  ht_indirect1 =  IE_sd(ypop = ht_ypop$ypop[2,], ygroup = yhat_group_sd2$yhat_group[,2,], 
                        hajofygroup = haj$haj[2,], boots = boots_est$boots, horvitzthompson = T, treatment = 2)[,,ngam]
  
  ht_oe = OE_sd(ypop = ht_ypop, ygroup = yhat_group_sd2$oe_yhat_group, 
                hajofygroup = haj$oe_haj, boots = boots_est$oehtboots, horvitzthompson = T)
  ht_oe_truevar = ht_oe$true_oe_var
  ht_oe = ht_oe$oe[,,ngam]
  
  ht_oe_mcvar = apply(boots_est$oe_est_ht, 1, var)
  haj_oe_mcvar = apply(boots_est$oe_est_ht,1, var)
  
  if(use.boot==T){
    
    return(list(unadj = po,
                haj = haj$haj,
                wtlist = haj$wts,
                direct = direct,
                ht_direct = ht_direct,
                indirect0 = indirect0,
                indirect1 = indirect1,
                ht_indirect0=ht_indirect0,
                ht_indirect1=ht_indirect1,
                ht_oe=ht_oe, #ht_oe_cov = ht_oe_cov,
                het_ie_truth = truth,
                oe = oe, oe_cov = oe_cov, ht_oe_truevar = ht_oe_truevar,
                ht_oe_mcvar = ht_oe_mcvar,
                haj_oe_mcvar = haj_oe_mcvar))
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

test = get_yhats_boot(clust_size, use.boot = T)

save.image(paste0('diffusion2/scenario', idx, concor, '.RSave'))


