#run this a single time and then use for all of the simulated datasets
library(gamlss)
library(rje)
library(tidyverse)
library(here)
print(here)
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

clust_avg_y = array(0, dim = c(n_clus, n_time, n_gamma, n_rep))

for(rep in 1:n_rep){ 
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
    
    # DRAW TREATMENT VECTOR USING DENSITY FROM LINEAR MODEL (NUMERATOR)
    sim_df$pi_itj = dBE(x = sim_df$A, mu = expit(sim_df$lin_pred))
    
    sim_df$new_A = rBE(n = nrow(sim_df), mu = expit(sim_df$lin_pred))
    
    sim_df = sim_df %>% group_by(neigh, time) %>% 
      mutate(new_W = mean(new_A), 
             new_U = mean(new_A[X1 == 1]), 
             new_V = mean(new_A[X2 == 1])) #should exclude unit j
    
    #OUTCOME MODEL
    sim_df$pot_out = beta_0 + beta_1*sim_df$new_A +  beta_2*sim_df$X1 + beta_3*sim_df$new_W + beta_4*sim_df$new_U + beta_5*sim_df$new_V 
    
    sim_df$pot_out = exp(sim_df$pot_out)
    
    neigh_results = sim_df %>% group_by(neigh, time) %>% summarise(Y = mean(pot_out))
    pop_results = sim_df %>% group_by(time) %>% summarise(Y = mean(pot_out))
    tot_results = sim_df %>% summarise(Y = mean(pot_out))
    
    #PUT IT ALL IN AN ARRAY
    for(nn in 1:n_clus){
      for(tt in 1:n_time){
        clust_avg_y[nn, tt, g, rep] = neigh_results %>% filter(neigh == nn, time == tt) %>% pull(Y)
    
      }
    }
  }
}

save(clust_avg_y, file = here('ch3truth', paste('ch3_truth', i, '.Rsave')))


#load('ch3_truth.Rsave')
#hist(sim_df$new_A)
