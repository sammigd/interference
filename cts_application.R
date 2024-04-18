library(gamlss)
library(stringr)
library(tidyverse)
library(rje)
library(here)
library(Interference, lib.loc = getwd())
library(magrittr)


here()
#load('/interference/cts_GroupIPW.R')
#load('/interference/cts_numerator.R')
#load('/interference/denominator_cts.R')
#load('/interference/cts_FromAlphaToRE.R')
options(dplyr.summarise.inform = FALSE)



load("/gpfs/gibbs/project/forastiere/sgd37/vax/vaccine_data_for_aim3.RSave")
df$neigh = as.numeric(factor(str_sub(df$census_geocode, 3, 5)))
df$time = as.numeric(df$month_ind)
dta = df

dta$A = dta$transformed_coverage
dta$Y = dta$case_count
dta$vax_rate = dta$vax2_count / dta$pop_count
dta$target = dta$age15_to_24
target_name = 'age15to24'


neigh_ind <- NULL
alpha_re_bound = 10


#get range of acceptable gammas
g_list = c()
for(cl in unique(dta$neigh)){
  #cl=1
  gamma_sub = dta %>% filter(neigh == cl)
  mod = gamlss(data = gamma_sub, vax_rate ~ target, family = 'BE')
  g_list = append(g_list, coef(mod)[2])
}
q1 = quantile(g_list,.10)
q2 = quantile(g_list,.90)


ngam = 11
gl = seq(from = q1, to = q2, length.out = 10)
gamma_list = rbind(rep(0, 10),
                   c(gl))
gamma_cts = cbind(gamma_list, rep(0, 2))

#####################################################
# GET LIST OF OBS FOR EACH CLUSTER, TIMEPOINT #######
#####################################################

if (is.null(neigh_ind)) {
  neigh_ind <- sapply(1 : max(dta$neigh), function(x) which(dta$neigh == x))
  if (typeof(neigh_ind) != 'list'){#if(const_size == T){
    ni2 = list()
    for(col in 1:ncol(neigh_ind)){
      ni2[[col]] = neigh_ind[,col]
    }
    neigh_ind = ni2
  }
}

#create a list of the observations that match to each time point
time_ind <- sapply(1:max(dta$time), function(x) which(dta$time == x))
time_ind_list = list()
for(col in 1:ncol(time_ind)){
  time_ind_list[[col]] = time_ind[,col]
}

n_neigh <- length(neigh_ind)
n_time <- length(time_ind_list)

######################################################
# BUILD DENOMINATOR MODEL FOR WEIGHTS ############################
######################################################
mod_a = gamlss(data = dta,
               A ~ re(fixed=~cov2lag2_prop + cov2lag2_prop_sq, random=~1|neigh), 
               family = BE)

mod_b = gamlss(data = dta, 
               A ~ re(fixed = ~cov2lag2_prop + cov2lag2_prop_sq + 
                               death_lag1 + case_lag1 + case_lag2 + 
                               cov2lag2_prop*death_lag1 + cov2lag2_prop_sq*death_lag1 + 
                               cov2lag2_prop_sq*case_lag1 + cov2lag2_prop*case_lag1 + 
                               cov2lag2_prop*case_lag2 + case_lag1*death_lag1 +
                               Median_Income + pct_nh_white + svi + 
                               pct_college + pct_male + age65_to_74 + #took out 75+ 
                               age55_to_64 + age45_to_54 + age35_to_44, 
                      random = ~1|neigh), 
               family = BE)

#calc msm weights
#marginalize over re values
f_int <- function(b, ii, mod){
  lin_pred <- logit(fitted(mod))[ii]
  prob_trt <- expit(lin_pred + b)
  return(prob_trt)
}

d_a <- c()
for(ii in 1:nrow(dta)){
  ans <- integrate(f_int, lower = -20, upper = 20, ii, mod_a)
  d_a = append(d_a, ans$value)
}

d_b <- c()
for(ii in 1:nrow(dta)){
  ans <- integrate(f_int, lower = -20, upper = 20, ii, mod_b)
  d_b = append(d_b, ans$value)
}

dta$pred_a = d_a#dBE(dta$A, mu = fitted(mod_a))
dta$pred_b = d_b#dBE(dta$A, mu = fitted(mod_b))

#TAKE PRODUCT OVER TIME POINTS
dta = dta %>%
  group_by(neigh) %>% 
  arrange(time) %>%
  mutate(ratio = pred_a / pred_b, 
         weights = cumprod(ratio)) #this is exactly the weights i would put into a msm

#TAKE PRODUCT OVER INDIVIDUALS IN A CLUSTER (BY TIME POINT)
denom_values = dta %>%
  group_by(neigh, time) %>%
  summarise(denom_weights = prod(weights), #this is the density of the entire observed trt vec for the cluster
            cluster_size = n_distinct(census_geocode)) #because sum over that same density for each unit? maybe dont need this line

dta = merge(dta, denom_values, by =c('neigh', 'time'))
n_neigh = length(unique(dta$neigh))

hajek_y_it = matrix(ncol = 0, nrow = n_neigh * n_time)
hajek_y_t = matrix(ncol = 0, nrow =  n_time)

hajek_w_it = matrix(ncol = 0, nrow = n_neigh * n_time)


for(gg in 1 : ncol(gamma_cts)){
  curr_gamma_numer = gamma_cts[,gg]
  ####################################################
  # GET INTERCEPTS FOR NUMERATOR MODEL ###############
  ####################################################
  alpha_by_cluster_by_time = dta %>% group_by(time, neigh) %>% summarise(alpha_vec = mean(vax_rate))
  dta$intercept_list = NA
  for(nn in 1:n_neigh){
    for (tt in 1:n_time){
      alpha_use = (alpha_by_cluster_by_time %>% filter(neigh == nn, time == tt) %>% pull(alpha_vec))[1]
      
      use_rows = intersect(time_ind_list[[tt]], neigh_ind[[nn]])
      
      #GET INTERCEPT FOR MONTHT CLUSTER NN
      lin_pred <- cbind(1, as.matrix(dta[use_rows,'target'])) %*% curr_gamma_numer
      
      re_alpha <- FromAlphaToRE(alpha = alpha_use, lin_pred = lin_pred,
                                alpha_re_bound = alpha_re_bound)
      
      dta$intercept_list[use_rows] <- re_alpha
    }
  }
  
  ######################################################
  # BUILD NUMERATOR MODEL FOR WEIGHTS ############################
  ######################################################
  design_mat = as.matrix(cbind(1, dta[, 'target']))
  dta$lin_pred = design_mat %*% curr_gamma_numer + dta$intercept_list
  
  #this just checks if intercepts is working
  dta %>% group_by(neigh, time) %>% summarise(mean(transformed_coverage),
                                              mean(expit(lin_pred)))
  
  #get counterfactual density
  dta$pi_itj = dBE(x = dta$vax_rate, mu = expit(dta$lin_pred)) #pi function in the numerator of wt 
  
  #' Now we have individual level denominators and numerators
  #' sum across the individuals  yhat_it (gamma)
  dta_w_hajek_it = dta %>%
    group_by(neigh, time) %>% arrange(neigh, time) %>% 
    #mutate(pi_itj = prod(pi_itj)) %>%
    summarise(sum_numeratorY = sum(pi_itj*Y), #sum across individuals
              sum_numerator = sum(pi_itj)) %>%
    merge(denom_values, by = c('neigh', 'time')) %>% #denom is already product the sum across individuals in the cluster
    mutate(hajek_num = sum_numeratorY / denom_weights,
           hajek_denom = sum_numerator / denom_weights,
           hajek_y = hajek_num / hajek_denom)
  hajek_y_it = cbind(hajek_y_it, dta_w_hajek_it$hajek_y)
  hajek_w_it = cbind(hajek_w_it, dta_w_hajek_it$hajek_denom) #shoudl be approx 1 or approx n?
  
  #sum accross the clusters to get a yhat_t(gamma)
  dta_w_hajek_t = dta_w_hajek_it %>% group_by(time) %>% arrange(time) %>%
    summarise(hajek_num = sum(hajek_num),
              hajek_denom = sum(hajek_denom)) %>%
    mutate(hajek_y = hajek_num / hajek_denom)
  hajek_y_t = cbind(hajek_y_t, dta_w_hajek_t$hajek_y)
}

#group avg hajek
haj_y_it_df = as.data.frame(hajek_y_it) %>%
  mutate( neigh = dta %>% group_by(neigh, time) %>% arrange(neigh, time) %>% summarise(x = 1) %>% pull(neigh),
          time = dta %>% group_by(neigh, time) %>% arrange(neigh, time) %>% summarise(x = 1) %>% pull(time))

#group w
haj_w_it_df = as.data.frame(hajek_w_it) %>%
  mutate( neigh = dta %>% group_by(neigh, time) %>% arrange(neigh, time) %>% summarise(x = 1) %>% pull(neigh),
          time = dta %>% group_by(neigh, time) %>% arrange(neigh, time) %>% summarise(x = 1) %>% pull(time))
#these weights are huge - they should be closer to n i think

#population average hajek estimate of Y
haj_y_t_df = as.data.frame(hajek_y_t) %>%
  mutate(time = dta %>% group_by(time) %>% arrange(time) %>% summarise(x = 1) %>% pull(time))


#now to plot....
obsy = dta %>% group_by(time, neigh) %>% summarise(y_it_obs = mean(case_count)) %>% group_by(time) %>% summarise(y_t_obs = mean(y_it_obs))
obsy = dta %>% group_by(time, neigh) %>% summarise(y_it_obs = mean(case_count)) 

plot_df <- haj_y_t_df %>%
  pivot_longer(cols = V1:V10, names_to = 'dirtygamma', values_to = 'Yhat') %>% 
  mutate(clean_gamma = factor(dirtygamma, levels = paste0('V', 1:10)),
         clean_gamma = recode(clean_gamma, !!!setNames(as.list(gamma_cts[2, 1:10]), paste0('V', 1:10))),
         clean_gamma = as.numeric(clean_gamma),
         oe_haj = Yhat - V11) %>%
  merge(obsy, by = 'time') %>%
  mutate(time = month(time, label = T),
         neigh = as.character(neigh))

ggplot(plot_df, aes(x = clean_gamma, y = Yhat)) + 
  geom_line(aes(y = y_it_obs, group = neigh, colour = neigh)) +
  geom_point() + 
  facet_wrap(~time, nrow = 1) + 
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = rgb(0,0,0,0)), 
        legend.position = 'none') + 
  ylab('Population Average Case Count') +
  xlab(expression(paste(gamma[''])))
ggsave(here('figures', 'ch3_figs', paste0('yhat_',target_name, '.png')), width = 6, height = 4)


ggplot(plot_df, aes(x = clean_gamma, y = oe_haj)) + 
  geom_point() + 
  facet_wrap(~time, nrow = 1) + 
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = rgb(0,0,0,0))) + 
  ylab('OE') +
  xlab(expression(paste(gamma[''])))
ggsave(here('figures', 'ch3_figs', paste0('oe_', target_name,'.png')), width = 6, height = 4)
