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
BootVar_sd <- function(dta, B = 500, alpha, ps = c('true', 'est'), cov_cols,
                       phi_hat_true = NULL, ps_info_est = NULL, verbose = TRUE,
                       ps_specs = NULL, trt_col = NULL, out_col = NULL,
                       return_everything = FALSE,
                       gamma_numer) {
  #B = 36
  num_gamma_cols <- ifelse(is.null(ncol(gamma_numer)), 1, ncol(gamma_numer))
  
  n_neigh <- max(dta$neigh)
  
  
  #re_var_positive <- rep(NA, B)
  library(pbmcapply)
  
  B = 36
  n_cores <- detectCores()
  output = pbmclapply(rep(1,n_cores), mc.cores = n_cores, function(i) {
    res <- replicate(ceiling(B/n_cores), inner_boot_funct(i, dta))
    res
  }) 
  
  #output = lapply(rep(1,10), function(i) {
  #  res <- replicate(1, inner_boot_funct(i, dta))
  #  res
  #}) 
  
  
  output = do.call(cbind, output)

  #restucture bootoutput
  hajboots = (abind(output['hajboots',], along = 3))#abind(output[[1]]['hajboots',], along = 3)
  ygroupboots = (abind(output['ygroup',], along = 3))#abind(output[[1]]['ygroup',], along = 3)
  chosenboots = (abind(output['chosen_clusters',], along = 2))#abind(output[[1]]['chosen_clusters',], along = 2)
  bootsboots = (abind(output['boots',], along = 3))#abind(output[[1]]['boots',], along = 3)
  oehajboots = (abind(output['oehajboots',], along = 2))#abind(output[[1]]['hajboots',], along = 3)
  
  
  return(list(hajboots = hajboots,
              ygroup = ygroupboots,
              chosen_clusters=chosenboots,
              boots = bootsboots,
              oehajboots = oehajboots)) 
} 
  
  
inner_boot_funct <- function(bb, dta){
  boot_dta <- GetBootSample(dta)
  chosen_clusters <- boot_dta$chosen_clusters
  
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
  
  ygroup <- ygb
  
  boots <- apply(ygb, c(2, 3), mean) 
  hajboots <- make_hajek(ygroup_boot, boot_dta)
  oehajboots <- hajboots$oe_haj
  hajboots <- hajboots$haj
  
  
  return(list(boots = boots, 
              ygroup = ygroup,
              chosen_clusters = chosen_clusters,
              hajboots = hajboots,
              oehajboots = oehajboots))
 
} 
  

