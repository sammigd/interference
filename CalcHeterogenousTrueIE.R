#create new function: get true vector of IE under heterogeneous interference 

get_het_ie <- function(dta, gamma_numer, cov_cols, interference = c('none', 'homog', 'hetero'),
                       beta_0, beta_1, beta_2, beta_3, beta_4, alpha, alpha_re_bound = 10){
  #beta_3  = 1
  #beta_2 = 0
  #interference = 'hetero'
  #beta_4 = 2
  clust_avg = array(0, dim = c(length(unique(dta$neigh)), ncol(gamma_numer), 2))
  clust_avg_oe = array(0, dim = c(length(unique(dta$neigh)), ncol(gamma_numer)))
  
  nn = 100
  for (sim in 1:nn){
    #sim = 1
    for (clus in unique(dta$neigh)){ #for each cluster
      #get covar and outcome vectors for cluster 
      #Xj = dta[dta$neigh == clus , cov_cols]
      #clus = 1
      
      Xj = dta[dta$neigh == clus , cov_cols]
    
      for (g in 1:ncol(gamma_numer)){
        #g = 1
        
        gamma_use = gamma_numer[,g]
        lin_pred <- cbind(1, as.matrix(Xj)) %*% gamma_use
        
        re_alpha <- FromAlphaToRE(alpha = alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        probs = expit(lin_pred + re_alpha)
        
        new_Aj = rbinom(length(probs),1,probs) 
        #redraw if no variation in Aj
        if(sum(new_Aj) == length(probs) | sum(new_Aj) == 0){new_Aj = rbinom(length(probs),1,probs)}
        
        
        new_Tj = sum(new_Aj / length(probs))
        #Tprimej = sum(new_Aj[Xj[,1] == 1]) / length(probs)
        Tprimej = sum(new_Aj[Xj[,1] == 1])/ sum(new_Aj)
        
        
        #outcome model
        pot_out = beta_0 + beta_1*new_Aj + 
          (as.matrix(Xj, ncol = 3) %*% c(beta_2, 0, 0)) + 
          #beta_2*Xj +
          beta_3*new_Tj + beta_4*Tprimej
          
        Yj_a0 = mean(pot_out[new_Aj == 0])
        Yj_a1 = mean(pot_out[new_Aj == 1])
        Yj = mean(pot_out)
          
        clust_avg[clus,g,1] = clust_avg[clus,g,1] + Yj_a0
        clust_avg[clus,g,2] = clust_avg[clus,g,2] + Yj_a1
        clust_avg_oe[clus,g] = clust_avg_oe[clus,g] + Yj
      }
    }
  }
  clust_avg = clust_avg / nn #aveg within each cluster
  clust_avg2 = apply(clust_avg, c(3,2), mean, na.rm = T)  #avg over all of the clusters
  
  clust_avg_oe = clust_avg_oe / nn
  clust_avg_oe2 =  apply(clust_avg_oe, 2, mean, na.rm = T)
  
  ie = array(NA, dim = c(2,length(gamma_numer[2,]),length(gamma_numer[2,])))
  oe = array(NA, dim = c(length(gamma_numer[2,]), length(gamma_numer[2,])))
  
  for (a1 in 1 : length(gamma_numer[2,])) {
    for (a2 in 1 : length(gamma_numer[2,])) {
      ie[, a1, a2] <- clust_avg2[,a1] - clust_avg2[,a2]
      oe[a1, a2] <- clust_avg_oe2[a1] - clust_avg_oe2[a2]
    }
  }

  
  return(list(clust_avg2, ie = ie, oe = oe,
              true_y_ie = clust_avg, true_y_oe = clust_avg_oe))
}

