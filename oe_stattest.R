oe_sigtest <- function(oe_output, oe_cov){ 
  #make sure you only input the 1-33 rows
  
  #testing
  #oe_output = oe
  #oe_cov = oe_cov[1:33, 1:33]

  #calculate test statistic
  oe_est = oe_output['est',1:33]
  oe_cov = oe_cov[1:33, 1:33]
  
  #get differences between each oe estimate
  diff_mat = as.matrix(dist(oe_est, method = 'manhattan')) #fix to not be just 1-3
  
  #get biggest diff
  max_diff_idx = which(diff_mat == max(diff_mat), arr.ind = TRUE)
  max_diff_val = diff_mat[max_diff_idx[1,1], max_diff_idx[1,2]] #this is the test statistic
  

  #simulate distribution of the test statistic
  test_stats = c()
  for (id in 1:1000){
    draw = MASS::mvrnorm(n = 1, mu = rep(0, nrow(oe_cov)), Sigma = oe_cov)
    
    #get differences between each oe estimate
    diff_mat_sim = as.matrix(dist(draw, method = 'manhattan'))
    
    #get biggest difference
    #max_diff_idx_sim = which(diff_mat_sim == max(diff_mat_sim), arr.ind = TRUE)
    #max_diff_val_sim = diff_mat[max_diff_idx_sim[1,1], max_diff_idx_sim[1,2]] #this is the test statistic
    max_diff_val_sim = max(diff_mat_sim)
    
    #add to list of test statistics
    test_stats = append(test_stats, max_diff_val_sim)
  }
  
  #this is quantile not ci
  ci = quantile(test_stats, c(0.025, 0.975))
  #make one sided
  ci_onesided = quantile(test_stats,.95)
 
  #return(list(accept = between(max_diff_val, ci[1], ci[2]), diff = max_diff_val, teststats = test_stats))
  return(list(accept = max_diff_val < ci_onesided, diff = max_diff_val, teststats = test_stats))
  
}
