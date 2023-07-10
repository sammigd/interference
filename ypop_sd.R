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
    } 
    
    
    #for de
    int2 = matrix(rep(0, (2*ngam)^2), nrow = (2*ngam))
    int2_oe = matrix(rep(0, (2*ngam)^2), nrow = (2*ngam))
    
    for (clus in 1:n_neigh){
      int1 = c(psi_i[clus,1,], psi_i[clus,2,]) #this is psi
      int1_oe = c(psi_oe[clus,], psi_oe[clus,])
      
      print(str(int1))
      print(str(int2))
      print(str((int1 %*% t(int1))))
      
      int2 = int2 + (int1 %*% t(int1)) #psi * t(psi) summed over clusters
      int2_oe = int2_oe + (int1_oe %*% t(int1_oe))
    }
      
    v = int2/n_neigh #this is expectation of score funct sd (sum over all ind than divide by n)
    a_de = c(wts[1,], wts[2,]) / n_neigh
    a_inv = solve(diag(a_de))
    sln = a_inv %*% v %*% t(a_inv) #now need to save for each gamma
    d_mat = cbind(diag(-1, nrow = ngam), diag(1, nrow = ngam))
    ypop_var_hj_de = (d_mat %*% sln %*% t(d_mat)) / n_neigh
    
    v_oe = int2_oe / n_neigh
    a_oe = c(ypop0$oe_wts, ypop0$oe_wts) / n_neigh
    a_oe_inv = solve(diag(a_oe))
    sln_oe = a_oe_inv %*% v_oe %*% t(a_oe_inv)
    
    #building ie0 d matrix
    m1 = matrix(rep(0, ngam^2), nrow =ngam)
    zero_matrix = matrix(rep(0, ngam^2), nrow =ngam)
    diag(m1) = -1
    m1[,ngam] <- m1[,ngam] + 1
    
    d_mat_ie0 = cbind(m1, zero_matrix)
    ypop_var_hj_ie0 = (d_mat_ie0 %*% sln %*% t(d_mat_ie0)) / n_neigh
    
    #building ie1 d matrix
    d_mat_ie1 = cbind(zero_matrix, m1)
    ypop_var_hj_ie1 = (d_mat_ie1 %*% sln %*% t(d_mat_ie1)) / n_neigh
    
    ypop_var_hj_ie = abind(ypop_var_hj_ie0, ypop_var_hj_ie1, along = 3)
    
    
    #for oe
    d_mat_oe = d_mat_ie0
    ypop_var_hj_oe = (d_mat_oe %*% sln_oe %*% t(d_mat_oe)) / n_neigh
    
    #tt = data.frame(x1 = oe$oe[2,,ngam], x2 = diag(ypop_var_hj_oe))
    
    if(FALSE){
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
    }
    
    return(list(ypop = ypop, ypop_var = ypop_var_hj, ypop_var_de = ypop_var_hj_de, ypop_var_ie = ypop_var_hj_ie, ypop_var_oe = ypop_var_hj_oe))
    
  }
}  
