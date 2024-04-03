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
  
  #pop_p_trt = sum(dta[,trt_col])/nrow(dta)
  
  # Specifyling neigh_ind will avoid re-running the following lines.
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
    #curr_gamma_numer = gamma_cts[,1]
    
    for (nn in 1 : n_neigh) {
      #nn = 1
      # Calculating the random effect that gives alpha.
      Xi <- dta[neigh_ind[[nn]], cov_cols]
      
      #not sure how to do RE for cts trt...
      if(F){
        lin_pred <- cbind(1, as.matrix(Xi)) %*% curr_gamma_numer
        re_alpha <- FromAlphaToRE(alpha = alpha, lin_pred = lin_pred,
                                alpha_re_bound = alpha_re_bound)}
      
      #some calcs for the overall effect
      Ai = dta[neigh_ind[[nn]], 'A']
      Ylagi = dta[neigh_ind[[nn]], c('case_lag1', 'death_lag1', 'death_lag2')]
      Alagi =  dta[neigh_ind[[nn]], c('cov2lag2_prop', 'cov2lag2_prop_sq')] 
      cl_df = df %>% filter(neigh == nn)
      oenum = cts_numerator(cl_df = cl_df, 
                            Ai_j = Ai, 
                            Xi_j = Xi, 
                            Alag = Alagi, 
                            Ylag = Ylagi,
                            gammas = curr_gamma_numer)
      
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
      
      for (curr_it in c(0, 1)) {
        #curr_it = 0
        bern_prob <- alpha ^ curr_it * (1 - alpha) ^ (1 - curr_it)
        prob_ind <- list(prob = 1)  
        y_curr <- 0
        wt_curr <- 0
        
        for (ind in neigh_ind[[nn]]) {
          #ind = 2
          if (dta$A[ind] == curr_it) {
            
            wh_others <- setdiff(neigh_ind[[nn]], ind)
            Ai_j <- dta$A[wh_others]
            Xi_j <- dta[wh_others, cov_cols]
            prob_ind <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j,
                                      gamma_numer = curr_gamma_numer,
                                      alpha = alpha, re_alpha = re_alpha)
            
            y_curr <- y_curr + dta$Y[ind] * prob_ind$prob
            wt_curr = wt_curr + (prob_ind$prob)
          }
        }
        
        denom <- calc_denominator_sd2(A = dta$A[neigh_ind[[nn]]], fix_phi = fix_phi, pop_p_trt = pop_p_trt)
        
        if(loud_denom == TRUE){print(paste0('denom is ', denom))}
        
        denom <- length(neigh_ind[[nn]]) * denom * bern_prob
        
        yhat_group[nn, curr_it + 1, gg] <- y_curr / denom
        wt_list[nn, curr_it + 1, gg] <- (wt_curr / denom) * length(neigh_ind[[nn]])
        
        #The denominator is the multiplied propensity scores of each individual in the cluster
        #that is, it is the propensity of the entire observed vector happening 
        
      }
    }
  }
  
  return(list(yhat_group = yhat_group, 
              wt_list = wt_list, 
              oe_yhat_group = oe_yhat_group, 
              oe_wt_list = oe_wt_list))
}
