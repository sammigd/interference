#aim 3

#' 
#' for each cluster we want weights that depend on gamma
#' 
#' want to calculate the treatment density from the observed data
#' this is the same as the denominator in the og version
#' 


cts_denom <- function(cl_df, 
                      Ai, 
                      Xi, 
                      Alag,
                      Ylag){
  #A: treatment vector
  #X: covariate matrix
  #svy_design_object: output of svydesign()
  
  #f_a = paste(A, '~', paste(Xi, collapse = '+'))
  #f_b = paste(A, '~', paste(c(Xi, Ylag), collapse = '+'))
  
  #little model
  mod_a = gamlss(Ai ~ Xi, family = BE)
  
  #big model
  df_b = cbind(Ai, Alag, Xi, Ylag)
  mod_b = gamlss(data = df_b, 
                 Ai ~ . , 
                 family = BE)

  #calc msm weights
  cl_df$pred_a = dBE(Ai, mu = fitted(mod_a))
  cl_df$pred_b = dBE(Ai, mu = fitted(mod_b))
  
  cl_df = cl_df %>%
    group_by(census_geocode) %>% 
    arrange(month_ind) %>%
    mutate(ratio = pred_a / pred_b, 
           weights = cumprod(ratio)) 
  #this gives a denominator weight for each observation - the propensity score as observed. 
  return(cl_df$weights)
}

