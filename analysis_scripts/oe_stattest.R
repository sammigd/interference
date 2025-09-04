oe_sigtest <- function(oe_output, oe_cov, gammas_idx = c(1:ncol(gamma_numer)), seed = 0){ 
  #make sure you only input the 1-33 rows
  if(seed!=0){set.seed(seed)}
  #testing
  #oe_output = oe
  #oe_cov = oe_cov[1:33, 1:33]

  #calculate test statistic
  oe_est = oe_output['est', gammas_idx]
  oe_cov = oe_cov[gammas_idx, gammas_idx]
  
  #get differences between each oe estimate
  diff_mat = as.matrix(dist(oe_est, method = 'manhattan')) #fix to not be just 1-3
  
  #get biggest diff
  max_diff_idx = which(diff_mat == max(diff_mat), arr.ind = TRUE)
  max_diff_val = diff_mat[max_diff_idx[1,1], max_diff_idx[1,2]] #this is the test statistic
  
  max_oe = max(oe_est)

  #simulate distribution of the test statistic
  test_stats = c()
  maxoe_test_stats = c() #for if the test stat is max oe instead of difference
  for (id in 1:1000){
    draw = MASS::mvrnorm(n = 1, mu = rep(0, nrow(oe_cov)), Sigma = oe_cov)
    
    #get differences between each oe estimate
    diff_mat_sim = as.matrix(dist(draw, method = 'manhattan'))
    
    #get biggest difference
    #max_diff_idx_sim = which(diff_mat_sim == max(diff_mat_sim), arr.ind = TRUE)
    #max_diff_val_sim = diff_mat[max_diff_idx_sim[1,1], max_diff_idx_sim[1,2]] #this is the test statistic
    max_diff_val_sim = max(diff_mat_sim)
   
    dist_max_oe = max(draw) #for if the test stat is just the max oe
    
    #add to list of test statistics
    test_stats = append(test_stats, max_diff_val_sim)
    maxoe_test_stats = append(maxoe_test_stats, dist_max_oe)
  }
  
  #this is quantile not ci
  ci = quantile(test_stats, c(0.025, 0.975))
  #make one sided
  ci_onesided = quantile(test_stats,.95)
  ci_max = quantile(maxoe_test_stats, 0.95)
  
  pvalue = sum(test_stats >= max_diff_val)/ length(test_stats)
 
  #return(list(accept = between(max_diff_val, ci[1], ci[2]), diff = max_diff_val, teststats = test_stats))
  return(list(accept = max_diff_val < ci_onesided, 
              diff = max_diff_val, 
              teststats = test_stats, 
              pvalue = pvalue,
              accept2 = max_oe < ci_max,
              maxoe = max_oe))
  
}
