#run this a single time and then use for all of the simulated datasets
library(gamlss)
library(rje)
library(tidyverse)
library(here)

source(here('interference', 'cts_functions.R'))
source(here('interference', 'helper_functs.R'))

options(dplyr.summarise.inform = FALSE)

args <- commandArgs(TRUE)
i = as.numeric(args[1])

# SET PARAMETERS
n_clus = 200
n_i <- 15
n_time <- 2
n_gamma <- 11
concor = 0.5
n_rep = 1

beta_0 = 3
beta_1 = 0
beta_2 = 0
beta_3 = 1
beta_4 = 1
beta_5 = 0

ngam = 11
gl = seq(from = -.8, to = .8, length.out = 10)
gamma_list = rbind(rep(0, 10),
                   c(gl))
gamma_cts = cbind(gamma_list, rep(0, 2))

clust_avg_y_sim = array(0, dim = c(n_clus, n_time, n_gamma, n_rep))
hajek_wts = array(0, dim = c(n_clus, n_time, n_gamma, n_rep))


for(rep in 1:n_rep){  
  #rep = 1
  print(rep)
  # GENERATE DATA USING PARAMETERS
  sim_df = make_simulated_dataset(n_clus, n_i, n_time, n_gamma, concor, 
                                  beta_0, beta_1, beta_2, beta_3, beta_4, beta_5, gamma_cts)
  
  #list of obs for each neigh
  neigh_ind <- sapply(1 : max(sim_df$neigh), function(x) which(sim_df$neigh == x))
  neigh_ind_list = list()
  for(col in 1:ncol(neigh_ind)){
    neigh_ind_list[[col]] = neigh_ind[,col]
  }
  
  #list of obs for each time
  time_ind <- sapply(1:max(sim_df$time), function(x) which(sim_df$time == x))
  time_ind_list = list()
  for(col in 1:ncol(time_ind)){
    time_ind_list[[col]] = time_ind[,col]
  }
  
  #DENOMINATOR MODEL FOR WEIGHTS
  mod_a = gamlss(data = sim_df,
                 A ~ re(fixed= ~ Alag + I(Alag^2), random=~1|neigh), 
                 family = BE)
  
  mod_b = gamlss(data = sim_df, 
                 A ~ re(fixed = ~Alag + I(Alag^2) + X1 + X2, #need to add in Ylag here
                        random = ~1|neigh), 
                 family = BE)
  
  f_int <- function(b, ii, mod){
    lin_pred <- logit(fitted(mod))[ii]
    prob_trt <- expit(lin_pred + b)
    return(prob_trt)
  }
  
  d_a <- c()
  for(ii in 1:nrow(sim_df)){
    ans <- integrate(f_int, lower = -20, upper = 20, ii, mod_a)
    d_a = append(d_a, ans$value)
  }
  
  d_b <- c()
  for(ii in 1:nrow(sim_df)){
    ans <- integrate(f_int, lower = -20, upper = 20, ii, mod_b)
    d_b = append(d_b, ans$value)
  }
  
  sim_df$pred_a = d_a#dBE(sim_df$A, mu = fitted(mod_a)) #d_a
  sim_df$pred_b = d_b#dBE(sim_df$A, mu = fitted(mod_b)) #d_b
  
  sim_df = sim_df %>%
    group_by(id) %>% #this needs to be indicidul level!!!!!
    arrange(time) %>%
    mutate(ratio = pred_a / pred_b, 
           weights = cumprod(ratio))
  
  #take product over time
  denom_values = sim_df %>%
    group_by(neigh, time) %>%
    summarise(denom_weights = prod(weights))
  
  sim_df = merge(sim_df, denom_values, by =c('neigh', 'time'))
  
  # GET TARGET INTERCEPT FOR NUMERATOR MODEL
  alpha_by_cluster_by_time = sim_df %>% group_by(time, neigh) %>% summarise(alpha_vec = mean(A))
  sim_df$intercept_list = NA
  
  for(g in 1:n_gamma){
    #g = 1
    curr_gamma_numer = gamma_cts[,g]
    for(nn in 1:n_clus){
      for(tt in 1:n_time){
        #nn = 1
        #tt = 1
        alpha_use = (alpha_by_cluster_by_time %>% filter(neigh == nn, time == tt) %>% pull(alpha_vec))[1]
        
        #PULL COVARIATES
        use_rows = intersect(time_ind_list[[tt]], neigh_ind_list[[nn]])
        #X <- sim_df[intersect(time_ind_list[[tt]], neigh_ind_list[[nn]]), c('X1', 'X2')] 
        
        #GET INTERCEPT 
        lin_pred <- cbind(1, as.matrix(sim_df[use_rows, c('X1')])) %*% curr_gamma_numer
        
        re_alpha <- FromAlphaToRE(alpha = alpha_use, lin_pred = lin_pred,
                                  alpha_re_bound = 10)
        
        sim_df$intercept_list[use_rows] <- re_alpha
      }
    }
    
    # FIT LINEAR MODEL FOR NUMERATOR
    design_mat = as.matrix(cbind(1, sim_df[, c('X1')]))
    sim_df$lin_pred = (design_mat %*% curr_gamma_numer)[,1] + sim_df$intercept_list
    
    sim_df %>% group_by(neigh, time) %>% summarise(mean(A),
                                                mean(expit(lin_pred)))
    
    #COUNTERFACTUAL NUMERATOR
    sim_df$pi_itj = dBE(x = sim_df$A, mu = expit(sim_df$lin_pred)) #pi function in the numerator of wt 
    
    #POTENTIAL OUTCOMES
    dta_w_hajek_it = sim_df %>%
      group_by(neigh, time) %>% arrange(neigh, time) %>% 
      summarise(sum_numeratorY = sum(pi_itj*Y), #sum across individuals
                sum_numerator = sum(pi_itj)) %>%
      merge(denom_values, by = c('neigh', 'time')) %>% 
      mutate(hajek_num = sum_numeratorY / denom_weights,
             hajek_denom = sum_numerator / denom_weights,
             hajek_y = hajek_num / hajek_denom)
    
    #PUT IT ALL IN AN ARRAY
    for(nn in 1:n_clus){
      for(tt in 1:n_time){
        clust_avg_y_sim[nn, tt, g, rep] = dta_w_hajek_it %>% filter(neigh == nn, time == tt) %>% pull(hajek_y)
        hajek_wts[nn, tt, g, rep] = dta_w_hajek_it %>% filter(neigh == nn, time == tt) %>% pull(hajek_denom)
        
        
      }
    }
  }
}

save(clust_avg_y_sim, file = here('ch3sim', paste0('ch3_sim', i, '.Rsave')))
save(hajek_wts, file = here('ch3sim_wts', paste0('ch3_simwt', i, '.Rsave')))



if(F){
  load(here('ch3_sim.Rsave')) #this has the truth
  load(here('ch3_truth.Rsave')) #this has the sim with the integration denom
  load(here('ch3_sim_withoutre.Rsave')) #this has the sim with the simple denom
  
  #get the mean sim estimates
  (apply(clust_avg_y_sim, c(3), mean)) #mean overall outcome
  (apply((clust_avg_y), c(3), mean))
  
  #simulation results
  test = data.frame(g = gamma_cts[2,-1],
                    y_haj = apply(clust_avg_y_sim, 3, mean),
                    y_haj_lb = apply(clust_avg_y_sim, 3, quantile, 0.025),
                    y_haj_ub = apply(clust_avg_y_sim, 3, quantile, 0.975),
                    truth = avg_truth[-1]) %>% #apply(clust_avg_y, 3, mean)) %>%
    mutate(bias = y_haj - truth)
  
  ggplot(test, aes(x = g, y = y_haj)) + 
    geom_point() + 
    geom_ribbon(aes(ymin = y_haj_lb, ymax = y_haj_ub, xmin = g, xmax = g), alpha = .2, fill = 'blue') +
    geom_line(aes(y = truth), colour = 'red')
  
  summary(sim_df$Y)
  
  #coverage and bias
  lb = apply(clust_avg_y_sim, 3, quantile, 0.025)
  ub = apply(clust_avg_y_sim, 3, quantile, 0.975)
  str(apply(clust_avg_y, c(3,4), mean))
  str(lb)
  
  coverage = array(NA, dim = c(10, n_rep))
  for(rep in 1:n_rep){
    true_y = apply(truth, c(3,4), mean)[,rep][-1]
    coverage[,rep] = between(true_y, lb, ub)
  }
  coverage = apply(coverage, 1, mean)
}

