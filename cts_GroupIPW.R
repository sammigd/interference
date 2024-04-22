#' Estimating the group average potential outcome
#' 
#' IPW estimator of the group average potential outcome.

GroupIPW_cts <- function(dta = df, 
                         cov_cols = c(13), 
                         gamma_numer = gamma_cts, 
                         alpha,
                         neigh_ind = NULL, 
                         trt_col = 26, 
                         out_col = 12, 
                         alpha_re_bound = 10, 
                         integral_bound = 10,
                         keep_re_alpha = FALSE,
                         verbose = TRUE,
                         loud_denom = FALSE, 
                         fix_phi = T, 
                         const_size = F, 
                         cts_trt = T) {
  
  #testing
  #alpha_re_bound = 15; verbose = F;gamma_numer = gamma_cts; loud_denom = F; fix_phi = T 
  
  integral_bound <- abs(integral_bound)
  alpha_re_bound <- abs(alpha_re_bound)
  dta <- as.data.frame(dta)
  
  # create a list of observations that match to each cluster
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
  
  yhat_group <- array(NA, dim = c(n_neigh, 2, ncol(gamma_numer)))
  dimnames(yhat_group) <- list(neigh = 1:n_neigh, trt = c(0, 1), gammas = 1:ncol(gamma_numer))
  
  oe_yhat_group <- array(NA, dim = c(n_neigh, ncol(gamma_numer)))
  dimnames(oe_yhat_group) <- list(neigh = 1:n_neigh, gammas = 1:ncol(gamma_numer))
  
  #create empty matrix for the hajek weights
  wt_list = array(0, dim = c(n_neigh, 2, ncol(gamma_numer)))
  dimnames(wt_list) <- list(neigh = 1:n_neigh, trt = c(0, 1), gammas = 1:ncol(gamma_numer)) #want one sum for each cluster
  # the sum of weights in each cluster
  
  #create empty matrix for the hajek weights
  oe_wt_list = array(0, dim = c(n_neigh, ncol(gamma_numer)))
  dimnames(oe_wt_list) <- list(neigh = 1:n_neigh, gammas = 1:ncol(gamma_numer))
  
  # Names of treatment and outcome column.
  if (!is.null(trt_col)) {names(dta)[trt_col] <- 'A'}
  if (!is.null(out_col)) {names(dta)[out_col] <- 'Y'}
  
  for (gg in 1 : ncol(gamma_numer)) {
    #gg = 1
    if (verbose){print(paste('gamma =', gamma_numer[-1,gg]))}
    
    curr_gamma_numer <- gamma_numer[,gg]

    #id target A for each cluster-timepoint
    alpha_by_cluster_by_time = dta %>% group_by(time, neigh) %>% summarise(alpha_vec = mean(A))
    
    #a single counterfactual for the whole population:
    Xi <- dta[, cov_cols]
    lin_pred = cbind(1, Xi) %*% curr_gamma_numer
    
    intercept_list = rep(NA, length(lin_pred))
    
    #a single denominator for the whole population
    oeden = cts_denom(cl_df = dta,
                      Ai = dta$A,
                      Xi = dta[,cov_cols],
                      Alag = dta[c('cov2lag2_prop', 'cov2lag2_prop_sq')],
                      Ylag = dta[c('case_lag1', 'death_lag1', 'death_lag2')],
                      neigh = dta$neigh)
    
    for (nn in 1 : n_neigh) {
      #nn = 1
      
      #FOR EACH CLUSTER AND EACH MONTH:
      #INSTEAD OF ALPHA, GET THE OBSERVED CLUSTER AVG TRT FOR EACH TIMEPOINT
      alpha_df = alpha_by_cluster_by_time %>% filter(neigh == nn)
      
      for (t in unique(dta$time)){
        t_rows = time_ind_list[[t]]
        
        alpha_use = (alpha_df %>% filter(time == t) %>% pull(alpha_vec))[1]
        
        #GET COVARIATES FOR MONTH T CLUSTER NN
        Xi <- dta[intersect(t_rows, neigh_ind[[nn]]), cov_cols] 
        
        #GET INTERCEPT FOR MONTHT CLUSTER NN
        lin_pred <- cbind(1, as.matrix(Xi)) %*% curr_gamma_numer
        
        re_alpha <- FromAlphaToRE(alpha = alpha_use, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        #ASSIGN THE RIGHT INTERCEPT TO EACH OBSERVATION
        intercept_list[intersect(time_ind_list[[t]], neigh_ind[[nn]])] <- re_alpha
      }
      
      #some calcs for the overall effect
      Ai = dta[neigh_ind[[nn]], 'A']
      Xi <- dta[neigh_ind[[nn]], cov_cols] 
      Ylagi = dta[neigh_ind[[nn]], c('case_lag1', 'death_lag1', 'death_lag2')]
      Alagi =  dta[neigh_ind[[nn]], c('cov2lag2_prop', 'cov2lag2_prop_sq')] 
      cl_df = df %>% filter(neigh == nn)
      oenum = cts_numerator(cl_df = cl_df, 
                            Ai_j = Ai, 
                            Xi_j = Xi, 
                            Alag = Alagi, 
                            Ylag = Ylagi,
                            gammas = curr_gamma_numer,
                            re_alpha = intercept_list[neigh_ind[[nn]]])
      
      ycurr = ((dta$Y[neigh_ind[[nn]]] * oenum)) #Yobs * counterfactual density
      
      oeden = cts_denom(cl_df = cl_df,
                        Ai = Ai,
                        Xi = Xi,
                        Alag = Alagi,
                        Ylag = Ylagi)
      
      oeden = length(neigh_ind[[nn]]) * oeden 
      
      wt_curr = (oenum)*length(neigh_ind[[nn]]) 
      #in the og version of the code oenum is a single product of trt probs
      #here do i want OE num to be a single numberwt
      
      oe_yhat_group[nn, gg] = sum(ycurr / oeden )
      oe_wt_list[nn, gg] <- wt_curr / oeden * length(neigh_ind[[nn]])  #(sum of all numerators / denom )* num of units in clus
      
      }
  }
  
  return(list(yhat_group = yhat_group, 
              wt_list = wt_list, 
              oe_yhat_group = oe_yhat_group, 
              oe_wt_list = oe_wt_list))
}
