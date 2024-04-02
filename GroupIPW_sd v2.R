#' Estimating the group average potential outcome
#' 
#' IPW estimator of the group average potential outcome.
#'
#' @param dta Data frame including treatment, outcome and covariates.
#' @param cov_cols The indices including the covariates of the ps model.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#' @param gamma_numer The coefficients of the ps model in the numerator.
#' If left NULL and estimand is 1, the coefficients in phi_hat will be used
#' instead.
#' @param alpha The values of alpha for which we want to estimate the group
#' average potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indices of
#' dta that are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'A' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#' @param integral_bound The number of standard deviations of the random effect
#' that will be used as the lower and upper limit.
#' @param keep_re_alpha Logical. If set to TRUE the "random" effect that makes
#' the average probability of treatment equal to alpha will be returned along
#' with the estimated group average potential outcome.
#' @param estimand Character, either '1' or '2.' If 1 is specified, then the
#' estimand with numerator depending on covariates is estimated. If estimand
#' is set equal to 2, the numerator considered is the product of Bernoulli.
#' @param verbose Whether printing of progress is wanted. Defaults to TRUE.
#' 
#' @export
GroupIPW_sd2 <- function(dta, cov_cols, gamma_numer = NULL, alpha,
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                     alpha_re_bound = 10, integral_bound = 10,
                     keep_re_alpha = FALSE, estimand = c('1', '2'),
                     verbose = TRUE, propensity_score = FALSE, loud_denom = FALSE, halfwt = F, 
                     fix_phi = T, const_size = F, cts_trt = F) {
  
  #testing
  #alpha_re_bound = 15; verbose = F;gamma_numer = gamma_list; loud_denom = F; fix_phi = T 
  
  estimand <- match.arg(estimand)
  integral_bound <- abs(integral_bound)
  alpha_re_bound <- abs(alpha_re_bound)
  dta <- as.data.frame(dta)
  
  pop_p_trt = sum(dta[,trt_col])/nrow(dta)
  
  # We only return the ksi's if we are estimating estimand 1.
  keep_re_alpha <- keep_re_alpha & (estimand == '1')
  
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
  
  #c reate empty matrix for the hajek weights
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
    
    for (nn in 1 : n_neigh) {
      #nn = 1
      # Calculating the random effect that gives alpha.
      Xi <- dta[neigh_ind[[nn]], cov_cols]
      lin_pred <- cbind(1, as.matrix(Xi)) %*% curr_gamma_numer
      re_alpha <- FromAlphaToRE(alpha = alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
      
      #some calcs for the overall effect
      Aj = dta[neigh_ind[[nn]], 'A']
      oenum = CalcNumerator(Ai_j = Aj, Xi_j = Xi,
                             gamma_numer = curr_gamma_numer,
                             alpha = alpha, re_alpha = re_alpha)
      
      ycurr = ((dta$Y[neigh_ind[[nn]]] * oenum$prob))
      
      oeden = calc_denominator_sd2(A = Aj, fix_phi = fix_phi, pop_p_trt = pop_p_trt)
      
      oeden = length(neigh_ind[[nn]]) * oeden 
      
      wt_curr = (oenum$prob)*length(neigh_ind[[nn]]) #sum of all the numerators in the cluster
      
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
