#' Bootstrap variance of potential outcomes.
#' 
#' Using re-sampling of clusters to acquire an estimate of the potential
#' outcome estimator variance.
#' 
#' @param dta The data frame including the observed data set.
#' @param B Number of bootstrap samples.
#' @param alpha The values of alpha where the potential outcomes are estimated.
#' @param ps Character. Whether the propensity score is known or is estimated.
#' Options include 'true', 'est'. Defaults to 'true'.
#' @param cov_cols Vector of column indices of the covariates used in the
#' propensity score model.
#' @param phi_hat_true List. Specify if the propensity score is known (ps set
#' to 'true'). Elements of the list are trt_coef and re_var including the
#' coefficients of the propensity score and random effect variance.
#' @param ps_info_est List of elements for acquiring estimates based on the
#' estimated propensity score. The list includes 1) glm_form: Element of 
#' formula class. The formula can be either for a fixed effects model or for a
#' model including random intercepts. 2) ps_with_re: An indicator of whether
#' the propensity score is a mixed model (set to TRUE) or not (set to FALSE).
#' 3) gamma_numer: The coefficients of the covariates in the counterfactual
#' treatment allocation model, and 4) use_control: Set to TRUE or FALSE if
#' you want or do not want additional elements in fitting the mixed model.
#' use_control does not have to be specified.
#' @param verbose Logical. Whether progress is printed. Defaults to TRUE.
#' @param trt_col If the treatment is not named 'A' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' @param return_everything Logical. Defaults to FALSE. If set to FALSE,
#' bootstrap estimates of the population average potential outcomes for all
#' values of alpha will be returned. If set to TRUE, additional information
#' will be returned including the chosen clusters at every bootstrap sample,
#' the group-specific estimated average potential outcome, and the random
#' effect variance of the cluster specific intercept.
#' 
#' @export
BootVar_sd <- function(dta, 
                       B = 500, alpha, 
                       gamma_numer = NULL,
                       #ps = c('true', 'est'), 
                       cov_cols,
                       phi_hat_true = NULL, 
                       ps_info_est = NULL, 
                       verbose = TRUE,
                       ps_specs = NULL, 
                       trt_col = NULL, 
                       out_col = NULL,
                       return_everything = FALSE,
                       phi_hat) {
  #B=10
  #horvitzthompson=F
  #ps <- match.arg(ps)
  n_neigh <- max(dta$neigh)
  
  chosen_clusters <- array(NA, dim = c(n_neigh, B))
  dimnames(chosen_clusters) <- list(neigh = 1 : n_neigh, sample = 1 : B)
  
  ygroup <- array(NA, dim = c(n_neigh, 2, ncol(gamma_numer), B))
  dimnames(ygroup) <- list(neigh = 1 : n_neigh, po = c('y0', 'y1'),
                           gamma = gamma_numer[2,], sample = 1 : B)
  
  boots <- array(NA, dim = c(2, ncol(gamma_numer), B))
  hajboots <- array(NA, dim = c(2, ncol(gamma_numer), B))
  
  dimnames(boots) <- dimnames(ygroup)[- 1]
  dimnames(hajboots) <- dimnames(boots)
  
  
  #re_var_positive <- rep(NA, B)
  
  for (bb in 1 : B) {
    
    if (verbose) {
      if (bb %% 10 == 0) {
        print(paste0('bootstrap sample ', bb))
      }
    }
    
    boot_dta <- GetBootSample(dta)
    chosen_clusters[, bb] <- boot_dta$chosen_clusters
    
    boot_dta <- boot_dta$boot_dta
    neigh_ind <- lapply(1 : max(boot_dta$neigh),
                        function(nn) which(boot_dta$neigh == nn))
    
    #assuming propensity score is true
      #re_var_positive[bb] <- (phi_hat_true[[2]] > 0)
      ygroup_boot <- GroupIPW_sd2(dta = boot_dta, cov_cols = cov_cols, phi_hat = phi_hat, 
                                  alpha = alpha, trt_col = trt_col, out_col = out_col,
                                  estimand = '1', alpha_re_bound = 20, verbose = F,
                                  gamma_numer = gamma_list,
                                  neigh_ind = neigh_ind,
                                  loud_denom = F, fix_phi = T, const_size = T)
      
      ygb = ygroup_boot$yhat_group
      
      ygroup[, , , bb] <- ygb
      
      boots[, , bb] <- apply(ygb, c(2, 3), mean) 
      hajboots[, , bb] <- make_hajek(ygroup_boot, boot_dta)$haj
    
      
  }
  
  if (return_everything) {
    return(list(boots = boots, ygroup = ygroup,
                chosen_clusters = chosen_clusters, hajboots = hajboots))
  }
  
  return(boots)
}
 
