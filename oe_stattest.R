oe_sigtest <- function(oe_output, oe_cov){ 
  #make sure you only input the 1-33 rows
  
  #testing
  #oe_output = oe
  #oe_cov = oe_cov[1:33, 1:33]

  #calculate test statistic
  oe_est = oe_output['est',]
  #diff_mat = as.matrix(dist(oe_est, method = 'manhattan')) #fix to not be just 1-3
  #max_diff_idx = which(diff_mat == max(diff_mat), arr.ind = TRUE)
  #max_diff_val = diff_mat[max_diff_idx[1,1], max_diff[1,2]] #this is the test statistic
  
  
  #restructure diff mat
  diff_tab = matrix(nrow = nrow(diff_mat)^2, ncol = 3)
  i = 1
  for(row in 1:nrow(diff_mat)){
    for (col in 1:ncol(diff_mat)){
      diff_tab[i,] <- c(diff_mat[row,col], row, col)
      i = i+1
    }
  }
  
  max_diff_val = max(abs(diff_tab[,1])) #abs not necceary bc include a-b and b-a
  
  #get new cov matrix corresponding to the differences
  new_cov = matrix(nrow = nrow(diff_tab), ncol = nrow(diff_tab))
  for (row1 in 1:nrow(diff_tab)){
    for (row2 in 1:nrow(diff_tab)){
      if(row1 == row2){ #calculate variances
        new_cov[row1,row2] <- oe_cov[diff_tab[row1,2], diff_tab[row1,2]] + 
                              oe_cov[diff_tab[row1,3], diff_tab[row1,3]] -
                              oe_cov[diff_tab[row1,2], diff_tab[row1,3]]
      }
      if(row1 != row2){ #calculate covariances
        a = diff_tab[row1,2]
        b = diff_tab[row1,3]
        c = diff_tab[row2,2]
        d = diff_tab[row2,3]
        new_cov[row1,row2] <- oe_cov[a,c] - oe_cov[a,d] - oe_cov[b,c] + oe_cov[b,d]
        
      }
    }
  }
  
  #isSymmetric(new_cov)
  #head(sort(eigen(new_cov)$values))

  #simulate distribution of the test statistic
  draws = MASS::mvrnorm(n = 1, mu = diff_tab[,1], Sigma = new_cov)
  ci = quantile(draws, c(0.025, 0.975))
  

  return(ci)

}
