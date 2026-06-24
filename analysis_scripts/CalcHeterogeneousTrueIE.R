#create new function: get true vector of IE under heterogeneous interference 
get_het_ie <- function(dta, gamma_numer, cov_cols, 
                       beta_0 = NULL, beta_1 = NULL, beta_2 = NULL, beta_3 = NULL, beta_4 = NULL, beta_5 = beta_5,
                       alpha, alpha_re_bound = 10, diffusion = F, diffusion_p = NULL, nn){

  clust_avg = array(0, dim = c(length(unique(dta$neigh)), ncol(gamma_numer), 2))
  clust_avg_oe = array(0, dim = c(length(unique(dta$neigh)), ncol(gamma_numer)))
  
  nn = nn 
  for (sim in 1:nn){
    print(paste0('start truth sim #', sim))
    for (clus in unique(dta$neigh)){ #for each cluster
      if(clus %% 50 ==0){print(paste0('cluster #', clus))}
      mini = dta %>% filter(neigh == clus)
      Xj = mini[,cov_cols]
    
      for (g in 1:ncol(gamma_numer)){
        gamma_use = gamma_numer[,g]
        lin_pred <- cbind(1, as.matrix(Xj)) %*% gamma_use
        
        re_alpha <- FromAlphaToRE(alpha = alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        probs = expit(lin_pred + re_alpha)
        new_Aj = rbinom(length(probs),1,probs) 
        
        #redraw if no variation in Aj
        if(sum(new_Aj) == length(probs) | sum(new_Aj) == 0){new_Aj = rbinom(length(probs),1,probs)}
        if(sum(new_Aj) == length(probs) | sum(new_Aj) == 0){new_Aj = rbinom(length(probs),1,probs)}
        
        if(diffusion == F){
          new_Tj = sum(new_Aj) / length(probs)
          Tprimej = sum(new_Aj*Xj[,1]) / sum(new_Aj)
          Tprimej2 = sum(new_Aj*Xj[,2]) / sum(new_Aj)
          
          #outcome model
          mini$pot_out = beta_0 + beta_1*new_Aj + 
            (as.matrix(Xj, ncol = 3) %*% c(beta_2, 0, 0)) + 
            beta_3*new_Tj + beta_4*Tprimej + beta_5*Tprimej2 + beta_6*new_Aj*ifelse(nrow(mini) == 15, 1, 0)
        }
        if(diffusion == T){
          mini$Aij = new_Aj
          #CALCULATE Y USING A
          mini$is_center_trted <- rep(mini %>% slice(1) %>% pull(Aij), each = 5)
          mini$Tij <- rep(mini %>% summarise(Tij = sum(Aij)) %>% pull(Tij), each = 5)

          mini = mini %>%
            mutate(n_untrt_alter = sum(Aij==0 & X1ij==0 & is_center_trted == 1))# %>%

          if(F){
            mini = mini %>%
              mutate(pot_out = case_when(Aij == 1 ~ 1, #treated node
                                         #untrt alter where center IS NOT treated
                                         Aij == 0 & X1ij == 0 & is_center_trted == 0 ~ 0,
                                         #untreated center - #max of a coin flip for each treated alter in the cluster
                                         Aij == 0 & X1ij == 1 ~ as.numeric(max(rbinom(Tij, 1, diffusion_p))), 
                                         #untrt alter where center IS treated
                                         Aij == 0 & X1ij == 0 & is_center_trted == 1 ~ as.numeric(rbinom(n_untrt_alter, 1, diffusion_p)))) #change the first 1 to number of rows
          }
          
          mini$pot_out[mini$Aij == 1] <- 1
          
          #untrt alter where center IS NOT treated
          mini$pot_out[mini$Aij == 0 & mini$X1ij == 0 & mini$is_center_trted == 0] <- 0
          
          #untreated center - #max of a coin flip for each treated alter in the cluster
          mini$pot_out[mini$Aij == 0 & mini$X1ij == 1] <- as.numeric(max(rbinom(mini$Tij[1], 1, diffusion_p)))
          
          #untrt alter where center IS treated
          mini$pot_out[mini$Aij == 0 & mini$X1ij == 0 & mini$is_center_trted == 1] <- as.numeric(rbinom(mini$n_untrt_alter, 1, diffusion_p))
        }

          
        Yj_a0 = mean(mini$pot_out[new_Aj == 0]) 
        Yj_a1 = mean(mini$pot_out[new_Aj == 1]) 
        Yj = mean(mini$pot_out)
          
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

