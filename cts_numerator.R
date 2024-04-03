cts_numerator <- function(cl_df, 
                          Ai_j, 
                          Xi_j,
                          Alag,
                          Ylag, 
                          gammas){
  
  design_mat_a = as.matrix(cbind(1, Ai_j))
  design_mat_b = as.matrix(cbind(1, Xi_j, Ai_j, Ylag))
  
  gammas_a = 0#gammas[1:2]
  gammas_b = gammas
  
  lin_pred_a = 0#design_mat_a %*% gammas_a #+ re_alpha_a
  lin_pred_b = design_mat_b %*% gammas_b #+ re_alpha_b
  
  cl_df$pred_a = dBE(x = Ai_j, mu = exp(lin_pred_a)-.01) #change to zero or 1 infl maybe???
  cl_df$pred_b = dBE(x = Ai_j, mu = exp(lin_pred_b)-.01)
  
  cl_df = cl_df %>%
    group_by(census_geocode) %>% 
    arrange(month_ind) %>%
    mutate(ratio = pred_a / pred_b, 
           weights = cumprod(ratio)) 
  return(cl_df$weights) 
}
