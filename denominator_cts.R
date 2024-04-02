#aim 3
library(gamlss)
#' 
#' for each cluster we want weights that depend on gamma
#' 
#' want to calculate the propensity score from the observed data
#' this is the same as the denominator in the og version
#' 
#' want to calcualte the propensity score from the counterfactual
#' this is what uses gammas
#' 
#' 
load("/gpfs/gibbs/project/forastiere/sgd37/vax/vaccine_data_for_aim3.RSave")

A = df$transformed_coverage
X = df[,13:20]
Ylag = df[,c('cov2lag2_prop', 'cov2lag2_prop_sq', 'death_lag1', 'death_lag2')]

A = 'transformed_coverage'
X = names(df)[13:20]
Ylag = c('cov2lag2_prop', 'cov2lag2_prop_sq', 'death_lag1', 'death_lag2')


cts_denom <- function(df, A, X, Ylag){
  #A: treatment vector
  #X: covariate matrix
  #svy_design_object: output of svydesign()
  f_a = paste(A, '~', paste(X, collapse = '+'))
  f_b = paste(A, '~', paste(c(X, Ylag), collapse = '+'))
  #little model
  mod_a = gamlss(data = df,
                 formula = as.formula(f_a), 
                 family = BE)
  
  mod_b = gamlss(data = df,
                 formula = as.formula(f_b), 
                 family = BE)
  
  #calc msm weights
  df$pred_a = dBE(df$transformed_coverage, mu = fitted(mod_a))
  df$pred_b = dBE(df$transformed_coverage, mu = fitted(mod_b))
  
  df = df %>%
    group_by(census_geocode) %>% 
    arrange(month_ind) %>%
    mutate(ratio = pred_a / pred_b, 
           weights = cumprod(ratio)) 
  #this gives a denominator weight for each observation - the propensity score as observed. 
  return(df$weights)
}
