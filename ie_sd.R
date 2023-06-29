#' Indirect effect estimates and asymptotic variance.
#' 
#' @param ygroup An matrix including the group average potential outcome
#' estimates where rows correspond to group, and columns to values of alpha.
#' #' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left
#' NULL for ps = 'true'.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @export
IE_sd <- function(ypop = NULL, ygroup, boots = NULL,
                  scores = NULL, alpha_level = 0.05, hajofygroup,
                  horvitzthompson = F, treatment = c(1,2)) {
  
  #testing
  #ygroup = yhat_group_sd2$yhat_group[,,1]
  #hajofygroup = haj$oe_haj
  #ps = 'true'
  
  #ypop = ypop 
  #ygroup = yhat_group_sd2$yhat_group[,1,]
  #hajofygroup = haj$haj[1,]
  #boots = boots_est$hajboots
  #treatment = 1
  
  gamma <- as.numeric(dimnames(ygroup)[[2]])
  ypop0 = ypop
  
  if(horvitzthompson == F){
    ypop <- hajofygroup
    ie_var = ypop0$ypop_var_ie[,,treatment]
  }else{
    ie_var = IEvar(ygroup = ygroup, ps = 'true', scores = scores)
  }
  
  names(ypop) <- gamma
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  
  dim_names <- c('est', 'var', 'LB', 'UB')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  ie <- array(NA, dim = c(length(dim_names), length(gamma), length(gamma)))
  dimnames(ie) <- list(stat = dim_names, gamma1 = gamma, gamma2 = gamma)
  
  for (a1 in 1 : length(gamma)) {
    for (a2 in 1 : length(gamma)) {
      ie[1, a1, a2] <- ypop[a1] - ypop[a2]
      #if(ypop[a1] - ypop[a2] > 1){print(paste(a1, a2))}
      #print(str(ie_var))
      ie[2, a1, a2] <- delta_method(ie_var[c(a1, a2), c(a1, a2)])
      ie_sd <- sqrt(ie[2, a1, a2])
      ie[c(3, 4), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
    }
  }
  ie_var_boots = NULL
  if (!is.null(boots)) {
    ie_var_boots <- array(NA, dim = c(length(gamma), length(gamma)))
    for (a1 in 1 : length(gamma)) {
      for (a2 in 1 : length(gamma)) {
        ie[5, a1, a2] <- var(boots[treatment, a1, ] - boots[treatment, a2, ])
        ie_sd <- sqrt(ie[5, a1, a2])
        ie[c(6, 7), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
        ie[c(8, 9), a1, a2] <- quantile(boots[1, a1, ] - boots[1, a2, ],
                                        probs = quants, na.rm = T)
      }
    }
  }
  
  return(list(ie = ie, ie_var = ie_var, ie_var_boots = ie_var_boots))
}
