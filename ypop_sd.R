#' Estimates and asymptotic variance of the population average potential
#' outcome for known or estimated propensity score model.
#' 
#' @param ygroup An array including the group average potential outcome
#' estimates where the dimensions correspond to group, individual treatment and
#' value of alpha.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left NULL
#' when ps is set to 'true'.
#' @param dta The data set including the variable neigh. Defaults to NULL. Can
#' be left NULL when the true propensity score is used.
#' @param use Whether the data with missing values will be used for the
#' estimation of the variance. Argument for cov() function. Defaults to
#' 'everything'.
#' 
#' @export
Ypop_sd <- function(ygroup, scores = NULL,
                    dta = NULL, use = 'everything', horvitzthompson = T) {
  
  #ygroup = yhat_group_sd2
  n_neigh <- dim(ygroup$yhat_group)[1]
  alpha <- as.numeric(dimnames(ygroup$yhat_group)[[3]])
  
  if(horvitzthompson == T){
    
    ypop <- apply(ygroup$yhat_group, c(2, 3), mean, na.rm = TRUE)
    oeypop = apply(ygroup$yhat_group, c(3), mean, na.rm = TRUE)

    ypop_var <- apply(ygroup$yhat_group, 3, function(x) cov(x, use = use))
    ypop_var <- array(ypop_var, dim = c(2, 2, ncol(gamma_numer)))
    # In order to get 1 / N, instead of 1 / (N - 1) in the variance estimates.
    ypop_var <- ypop_var * (n_neigh - 1) / n_neigh
    ypop_var <- ypop_var / n_neigh  # Since we have n_neigh clusters.
    dimnames(ypop_var) <- list(a = c(0, 1), a = c(0, 1), gamma = gamma_numer[2,])
    
    size = dta %>% group_by(neigh) %>% summarise(cluster_pop = n())
    psi_oe_ht = ygroup$oe_yhat_group
    
    for (clus in 1:n_neigh){
      c_size = unlist(size[clus, 'cluster_pop'])
      psi_oe_ht[clus,] = (ygroup$oe_yhat_group[clus,] * c_size) - (oeypop * c_size)
    }
    
    numgam = dim(ygroup$yhat_group)[3]
    A_inv = solve(diag(rep(-1, numgam)))
    
    int1 = apply(psi_oe_ht, 2, sum) 
    int2 = (int1 %*% t(int1))
    V =  int2 / nrow(size)#apply(psi_oe_ht^2, 2, mean)
    
    oe_ht_var = A_inv %*% V %*% t(A_inv)
    
    return(list(ypop = ypop, ypop_var = ypop_var, oe_ht_var = oe_ht_var))
    
  } else{  #if hajek
    
    #the estimates
    ypop0 <- make_hajek(ygroup, dta)
    wts <- ypop0$wts #summed over ind and clus
    ypop = ypop0$haj

    #the variance calc
    size = dta %>% group_by(neigh) %>% summarise(cluster_pop = n())
    a =  wts / n_neigh #wts here are summed across ind and clus
    
    #hajek at the cluster level
    haj_i = ygroup$yhat_group / ygroup$wt_list #apply(ygroup$yhat_group, c(2,3), "*", size$cluster_pop) / ygroup$wt_list

    psi_i = ygroup$yhat_group
    psi_oe = ygroup$oe_yhat_group

    for (clus in 1:n_neigh){
      psi_i[clus,,] = (ygroup$yhat_group[clus,,] * unlist(size[clus, 'cluster_pop'])) - (ypop * ygroup$wt_list[clus,,])#(haj_i[clus,,] * ygroup$wt_list[clus,,] * unlist(size[clus, 'cluster_pop'])) #this yhat should be the hajek estimates (not summed over cluster)
      psi_oe[clus,] = (ygroup$oe_yhat_group[clus,] * unlist(size[clus, 'cluster_pop'])) - (ypop0$oe_haj * ygroup$oe_wt_list[clus,])
    }
    
    #for ypop var
    ypop_var_hj <- array(NA, dim = c(2, 2, ncol(gamma_numer)))
    for (g in 1:dim(ygroup$yhat_group)[3]){
      int1 = psi_i[,,g] #the psi for a specific gamma
      int2 = matrix(c(0,0,0,0), nrow = 2)
      
      for (clus in 1:n_neigh){
        int2 = int2 + (int1[clus,] %*% t(int1[clus,]))
      }
      
      v = int2/n_neigh
      a_inv = solve(diag(a[,g]))
      
      sln = a_inv %*% v %*% t(a_inv) #now need to save for each gamma
      ypop_var_hj[,,g] = sln / n_neigh
    } #this and de version seem equivalent...
    
    
    #for de
    ypop_var_hj_de <- array(NA, dim = c(ncol(gamma_numer), ncol(gamma_numer)))
    int2 = matrix(rep(0, 128^2), nrow = 128)
      
    for (clus in 1:n_neigh){
      int1 = c(psi_i[clus,1,], psi_i[clus,2,]) #this is psi
      int2 = int2 + (int1 %*% t(int1)) #psi * t(psi) summed over clusters
    }
      
    v = int2/n_neigh
    a_de = c(wts[1,], wts[2,]) / n_neigh
    a_inv = solve(diag(a_de))
      
    sln = a_inv %*% v %*% t(a_inv) #now need to save for each gamma
    ypop_var_hj_de[,] = sln[ncol(gamma_numer) + 1 : ncol(gamma_numer*2), 1: ncol(gamma_numer)] / n_neigh
  
    #for the indirect effect hajek variance
    ypop_var_hj_gammalevel = array(NA, dim = c(ncol(gamma_numer), ncol(gamma_numer), 2))
    for (trt in 1:2){
      int1 = psi_i[,trt,]
      int2 = matrix(rep(0, ncol(gamma_numer)^2), nrow = ncol(gamma_numer))
      
      for (clus in 1:n_neigh){
        int2 = int2 + (int1[clus,] %*% t(int1[clus,]))
      }
      
      v = int2/n_neigh
      a_inv = solve(diag(a[trt,]))
      
      sln = a_inv %*% v %*% t(a_inv) 
      ypop_var_hj_gammalevel[,,trt] = sln / n_neigh
    }
    
    #for the overall effect hajek variance
    ypop_var_hj_oe = array(NA, dim = c(ncol(gamma_numer), ncol(gamma_numer)))
    int1 = psi_oe
    int2 = matrix(rep(0, ncol(gamma_numer)^2), nrow = ncol(gamma_numer))

    for(clus in 1:n_neigh){
      int2 = int2 + (int1[clus,] %*% t(int1[clus,]))
    }

    v = int2/n_neigh

    a = ypop0$oe_wts / n_neigh
    a_inv = solve(diag(a))

    sln = a_inv %*% v %*% t(a_inv)
    ypop_var_hj_oe = sln / n_neigh
    
    return(list(ypop = ypop, ypop_var = ypop_var_hj, ypop_var_hj_de = ypop_var_hj_de, ypop_var_ie = ypop_var_hj_gammalevel, ypop_var_oe = ypop_var_hj_oe))
    
  }
}  
