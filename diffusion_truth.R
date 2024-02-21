####################################################

#all possible trt vectors for cluster of 5
#trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1)) #each row is a possible trt vector

get_true_diffusion <- function(gamma_numer,
                               x2_prev = 0.2,
                               kappa = 0.25,
                               trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1)),
                               clust_size = 5,
                               p_diffusion = 0.5,
                               p_trt = 0.25){
  #SETUP RESULTS TABLES
  #P(A_i |X_i)
  res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(trt_vecs), ncol(gamma_numer)))
  dimnames(res_mat) <- list(trt_vecs = 1:ncol(trt_vecs), x_vec = 1:ncol(trt_vecs), gammas = 1:ncol(gamma_numer))
  
  #TEMPS 2/7
  x1_res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer) - 33))
  x1_ind_res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer) - 33, clust_size))
  x1_adj_res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer) - 33, clust_size))
  adj_Ybar = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer)-33,3))
  

  #E[Y_i(a)]
  Y_mat = array(NA, dim = c(ncol(trt_vecs), 3))
  dimnames(Y_mat) <- list(trt_vecs = 1:ncol(trt_vecs), y = c('Y', 'Y0', 'Y1'))
  
  #P(A_ij |X_ij)
  ind_res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(trt_vecs), ncol(gamma_numer), clust_size))
  dimnames(ind_res_mat) <- list(trt_vecs = 1:ncol(trt_vecs), x_vec = 1:ncol(trt_vecs), gammas = 1:ncol(gamma_numer), j = 1:clust_size)
  
  #E[Y_ij]
  ind_Y_mat = array(NA, dim = c(ncol(trt_vecs), clust_size))
  dimnames(ind_Y_mat) <- list(trt_vecs = 1:ncol(trt_vecs), j = 1:clust_size)
  
  #FIXED X1 FOR EVERYONE: ALTER/CENTER INDICATOR
  x1 = rep(0, clust_size); x1[1] = 1
  print(x1)
  
  print('*')
  
  #LOOP OVER POSSIBLE TREATMENT VECTORS
  for (a in 1:ncol(trt_vecs)){
    #a = 5
    aa = trt_vecs[,a]
      
    #INDIVIDUAL AVG EXPECTED OUTCOMES E[Y_ij]
    Yij = ifelse(aa == 1, 1, NA)
    Yij[1] = ifelse(aa[1] == 0, 1 - (1-p_diffusion)^sum(aa[2:clust_size]), Yij[1])
    Yij[2:clust_size] = ifelse(aa[1] == 1 & aa[2:clust_size] == 0, p_diffusion, Yij[2:clust_size])
    Yij[2:clust_size] = ifelse(aa[1] == 0 &  aa[2:clust_size] == 0, 0, Yij[2:clust_size])
    
    ind_Y_mat[a,] <- Yij
    
    #CLUSTER AVG EXPECTED OUTCOMES E[Y_i]
    Y_bar1 = 1
    if(aa[1] == 1){ #if the center is treated
      Y_bar = (sum(aa) + (clust_size - sum(aa))*p_diffusion) / clust_size
      Y_bar0 =( 4*p_diffusion + (1- (1-p_diffusion)^sum(aa[-1]))) / clust_size
      
    } 
    if(aa[1] == 0){ # if center is not treated
      Y_bar = (sum(aa) + (1- (1-p_diffusion)^sum(aa))) / clust_size
      Y_bar0 =( 0*p_diffusion + (1- (1-p_diffusion)^sum(aa[-1]))) / clust_size
    }
    
    Y_mat[a,] = c(Y_bar, Y_bar0, Y_bar1)
    
    #LOOP OVER GAMMAS WHERE INTERVENTION IS ON X1
    for (g in 1:(ncol(gamma_numer)-33)){
      #g = 14
      curr_gamma = gamma_numer[,g]
      
      lin_pred <- cbind(1, as.matrix(cbind(x1, xx))) %*% curr_gamma
      re_alpha <- FromAlphaToRE(alpha = p_trt, lin_pred = lin_pred,
                                alpha_re_bound = 15)
      
      ind_p_a_x = expit(re_alpha + lin_pred)^aa * (1-expit(re_alpha + lin_pred))^(1-aa)
      
      p_a_x = CalcNumerator(Ai_j = aa, 
                            Xi_j = cbind(x1, xx),
                            gamma_numer = curr_gamma,
                            alpha = p_trt, 
                            re_alpha = re_alpha, 
                            include_alpha = FALSE)$prob #equal to prod(ind_p_a_x)
      
      x1_res_mat[a,g] <- p_a_x
      x1_ind_res_mat[a,g,] <- ind_p_a_x
      #x1_adj_res_mat[a,g,] <- p_a_x / ind_p_a_x #changing gamma changes P(A|X very little)
      
      #CALCULATE ADJUSTED E[Yi-j] BY DIVIDING OUT P(Aij | Xij)
      adj_Yij = ind_Y_mat / x1_ind_res_mat[a,g,]
      use_Yij = mean(adj_Yij[a,])
      adj_Yi0  = ifelse(sum(aa) == clust_size, NA, mean(adj_Yij[a,][aa == 0])) #problem
      adj_Yi1  = ifelse(sum(aa) == 0, NA, mean(adj_Yij[a,][aa == 1])) #problem
      adj_Ybar[a,g,] <- c(use_Yij, adj_Yi0, adj_Yi1)
      }
  }
  Exp_Y = apply(x1_res_mat*adj_Ybar[,,1],2, sum, na.rm = T) #summing over all trt vecs
  Exp_Y0 = apply(x1_res_mat*adj_Ybar[,,2],2, sum, na.rm = T) 
  Exp_Y1 = apply(x1_res_mat*adj_Ybar[,,3],2, sum, na.rm = T) 
  
  Exp_Y1 - Exp_Y0
  
  
    
    #MULTIPLY P(A_i-j | X_i-j) by E[Yij]
    
    #LOOP OVER POSSIBLE X2 VECTORS
  print('**')
  print(ncol(trt_vecs))
    for (x in 1:ncol(trt_vecs)){
      print(x)
      #x = 1
      xx = trt_vecs[,x]
      print(xx)
      
      #CALCULATE P(X2|X1) - ONLY NEEDED FOR X2 INTERVENTION
      #p(x2|x1)
      n_central_x2_1 = sum(xx[1] == 1)
      n_central_x2_0 = sum(xx[1] == 0)
      n_alter_x2_1 = sum(xx[-1] == 1)
      n_alter_x2_0 = sum(xx[-1] == 0) #these n's must sum to five
     
      #n_central_x2_1 + n_central_x2_0 + n_alter_x2_1 + n_alter_x2_0
      
      p11 = (x2_prev + kappa)
      p01 = 1-p11
      p10 = x2_prev + (kappa / (length(aa) -1))
      p00 = 1 - p10
      
      p_x2_x1 = (p11^n_central_x2_1) * (p01^n_central_x2_0) * (p10^n_alter_x2_1) * (p00^n_alter_x2_0) 
      
      ind_p_x2_x1 = vector(length = 5)
      ind_p_x2_x1[1] = ifelse(xx[1] == 0, p01, p11)
      ind_p_x2_x1[2:clust_size] = ifelse(xx[2:clust_size] == 0, p00, p10)
      
      #LOOP OVER POSSIBLE GAMMAS
      for (g in 1:ncol(gamma_numer)){
        #g = 1
        curr_gamma = gamma_numer[,g]
        
        lin_pred <- cbind(1, as.matrix(cbind(x1, xx))) %*% curr_gamma
        re_alpha <- FromAlphaToRE(alpha = p_trt, lin_pred = lin_pred,
                                  alpha_re_bound = 15)
        
        ind_p_a_x = expit(re_alpha + lin_pred)^aa * (1-expit(re_alpha + lin_pred))^(1-aa)
        
        p_a_x = CalcNumerator(Ai_j = aa, 
                              Xi_j = cbind(x1, xx),
                              gamma_numer = curr_gamma,
                              alpha = p_trt, 
                              re_alpha = re_alpha, 
                              include_alpha = FALSE)$prob
        
        if(curr_gamma[3] == 0){
          res_mat[a,x,g] = p_a_x
          ind_res_mat[a,x,g,] <- ind_p_a_x
          
          }
        if(curr_gamma[3] != 0){
          res_mat[a,x,g] = p_a_x * p_x2_x1 #this is P(Aij | X2)
          ind_res_mat[a,x,g,] <- ind_p_a_x * ind_p_x2_x1 
          }
      }
    }
  
  ind_p_a_x1 = apply(ind_res_mat, c(1,3,4), sum, na.rm =T) #this gives a number for each individual in each treatment vector for each gamma
  #ind_p_a_x1 = ind_res_mat[,1,,]
  # ind_p_a_x1 = data.frame(fix_Aij = c(rep(0, 4), rep(1,4)),
  #                         fix_X1ij = rep(c(0,1), 4),
  #                         fix_X2ij = rep(c(0,0,1,1), 2)) %>%
  #   mutate(p_x2_x1 = (p11^(fix_X2ij == 1 & fix_X1ij == 1)) * 
  #                    (p01^(fix_X2ij == 0 & fix_X1ij == 1)) * 
  #                    (p10^(fix_X2ij == 1 & fix_X1ij == 0)) * 
  #                    (p00^(fix_X2ij == 0 & fix_X1ij == 0)) )
  # 
  
  
  #at clus level:
  #IF X2 INTERVENTION, NEED TO DEAL WITH X2
  p_a_x1 = apply(res_mat, c(1,3), sum) #P(A | X1) 
  
  #IF X1 INTERVENTION, DO NOT NEED TO TAKE WEIGHTED SUM! ALREADY HAVE P(A|X1) IN RES_MAT
  res_mat_x1 = res_mat[,1,34:67] #this is p(A_i | X1_i)
  ind_res_mat_x1 = ind_res_mat[,1,34:67,] #choose any x2
  
  #create the df for Ai-j | Xi-j
  p_aisubj_xisubj = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer)-33, clust_size))
  names(p_aisubj_xisubj) = names(ind_res_mat)[-2]
  
  for (j in 1:clust_size){
    p_aisubj_xisubj[,,j] = res_mat_x1 / ind_res_mat_x1[,,j]
  }
  
  x1_final = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer)-33, 3))
  for(t_vec in 1:ncol(trt_vecs)){
    Yi = Y_mat[t_vec,1]
    Y0i = Y_mat[t_vec,2]
    Y1i = Y_mat[t_vec,3] #possuble prob: Y1i = 1 even when no one in cluster is treated (t_Vec = 1)
    
    x1_final[t_vec,,1] = apply(p_aisubj_xisubj[t_vec,,] * Yi, 1, mean)  #have entry for each gamma 
    x1_final[t_vec,,2] = apply(p_aisubj_xisubj[t_vec,,] * Y0i, 1, mean)
    x1_final[t_vec,,3] = apply(p_aisubj_xisubj[t_vec,,] * Y1i, 1, mean)
  }
  x1_final = apply(x1_final, c(2,3), sum)

  #now just want to sum across all the treatment vectors to get a single number for each gamma
  
  #get yi's / ind_Y_mat
  adj_Ybar = array(NA, dim = c(ncol(trt_vecs), ncol(gamma_numer),3))
  for (g1 in 1 : length(gamma_numer[2,])) {
    adj_Yij = ind_Y_mat / ind_p_a_x1[,g1,]
    for (a in 1:ncol(trt_vecs)){
      aa = trt_vecs[,a]
      use_Yij = mean(adj_Yij[a,])
      adj_Yi1  = ifelse(sum(aa) == 0, NA, mean(use_Yij[aa == 1]))
      adj_Yi0  = ifelse(sum(aa) == clust_size, NA, mean(use_Yij[aa == 0]))
      adj_Ybar[a,g1,] <- c(use_Yij, adj_Yi0, adj_Yi1)
    }
  }
  
  #multiply adjusted Ybars by P(Ai = Xi)
  wtd_y_components = p_a_x1 * adj_Ybar[,,1]
  wtd_y = apply(wtd_y_components, 2, sum)
  
  wtd_y0_components = p_a_x1 * adj_Ybar[,,2]
  wtd_y0 = apply(wtd_y0_components, 2, sum, na.rm =T)
  
  wtd_y1_components = p_a_x1 * adj_Ybar[,,3]
  wtd_y1 = apply(wtd_y1_components, 2, sum, na.rm = T)
  
  
  ie = array(NA, dim = c(2, ncol(gamma_numer), ncol(gamma_numer)))
  oe = array(NA, dim = c(ncol(gamma_numer), ncol(gamma_numer)))
  
  for (g1 in 1 : length(gamma_numer[2,])) {
    for (g2 in 1 : length(gamma_numer[2,])) {
      ie[, g1, g2] <- c(wtd_y0[g1] - wtd_y0[g2], wtd_y1[g1] - wtd_y1[g2])
      oe[g1, g2] <- wtd_y[g1] - wtd_y[g2]
    }
  }
  
  de = wtd_y1 - wtd_y0
  
  
  return(list(de = de, ie = ie, oe = oe, y0 = wtd_y0, y1 = wtd_y1))
}

trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1))
results = get_true_diffusion(gamma_numer = gamma_numer, x2_prev = 0.2, kappa = 0.25, trt_vecs = trt_vecs)

