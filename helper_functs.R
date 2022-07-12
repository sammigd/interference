FromAlphaToRE <- function(alpha, lin_pred, alpha_re_bound = 10) {
  
  alpha_re_bound <- abs(alpha_re_bound)
  
  r <- optimise(f = AlphaToBi, lower = - 1 * alpha_re_bound,
                upper = alpha_re_bound, alpha = alpha, lin_pred = lin_pred)
  r <- r$minimum
  if (alpha_re_bound - abs(r) < 0.1) {
    warning(paste0('bi = ', r, ', alpha_re_bound = ', alpha_re_bound))
  }
  return(r)
}

AlphaToBi <- function(b, alpha, lin_pred) {
  exp_lin_pred <- exp(lin_pred)
  r <- abs(mean(exp_lin_pred / (exp_lin_pred + exp(- b))) - alpha)
  return(r)
}

CalcNumerator <- function(Ai_j, Xi_j, gamma_numer, alpha, re_alpha) {
  
  gamma_numer <- matrix(gamma_numer, nrow = length(gamma_numer), ncol = 1)
  
  lin_pred <- cbind(1, as.matrix(Xi_j)) %*% gamma_numer
  lin_pred <- lin_pred + re_alpha
  probs <- expit(lin_pred)
  
  r <- (probs / alpha) ^ Ai_j * ((1 - probs) / (1 - alpha)) ^ (1 - Ai_j)
  return(list(prob = prod(r), re_alpha = re_alpha))
}


GetBootSample <- function(dta) {
  
  num_clus <- max(dta$neigh)
  boot_clusters <- sample(1 : num_clus, num_clus, replace = TRUE)
  
  # Binding data without accidentally merging repeated clusters.
  boot_dta <- NULL
  for (nn in 1 : num_clus) {
    D <- subset(dta, neigh == boot_clusters[nn])
    D$neigh <- nn
    boot_dta <- rbind(boot_dta, D)
  }
  
  return(list(boot_dta = boot_dta, chosen_clusters = boot_clusters))
}
