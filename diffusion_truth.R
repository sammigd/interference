####################################################

#all possible trt vectors for cluster of 5
#trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1)) #each row is a possible trt vector

get_true_diffusion <- function(gamma_numer,
                               x2_prev = 0.2,
                               kappa = 0.25,
                               trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1)),
                               clust_size = 5,
                               p_diffusion = 0.5){
  #setup results tables
  res_mat = array(NA, dim = c(ncol(trt_vecs), ncol(trt_vecs), ncol(gamma_numer)))
  dimnames(res_mat) <- list(trt_vecs = 1:ncol(trt_vecs), x_vec = 1:ncol(trt_vecs), gammas = 1:ncol(gamma_numer))
  
  Y_mat = array(NA, dim = c(ncol(trt_vecs), 3))
  dimnames(Y_mat) <- list(trt_vecs = 1:ncol(trt_vecs), y = c('Y', 'Y0', 'Y1'))
  
  x1 = rep(0, clust_size); x1[1] = 1
  
  for (a in 1:ncol(trt_vecs)){
    #aa = 1
    aa = trt_vecs[,a]
    
    #Y_i
    Y_bar1 = 1
    if(aa[1] == 1){ #if the center is treated
      Y_bar = (sum(aa) + (clust_size - sum(aa))*p_diffusion) / clust_size
      Y_bar0 =( 4*p_diffusion + (1- (1-p_diffusion)^sum(aa[-1]))) / clust_size
      
      #Y_bar0 = p_diffusion
    } 
    if(aa[1] == 0){ # if center is not treated
      Y_bar = (sum(aa) + (1- (1-p_diffusion)^sum(aa))) / clust_size
      #Y_bar0 = (1 - (1 - p_diffusion)^sum(aa)) / (clust_size) #- sum(aa))
      Y_bar0 =( 0*p_diffusion + (1- (1-p_diffusion)^sum(aa[-1]))) / clust_size
      
    }
    
    #Y_bar0 =( 4*p_diffusion + (1- (1-p_diffusion)^sum(aa[-1]))) / clust_size
    
    Y_mat[a,] = c(Y_bar, Y_bar0, Y_bar1)
    
    for (x in 1:ncol(trt_vecs)){
      #xx = 1
      xx = trt_vecs[,x]
      
      #p(x2|x1)
      n_central_x2_1 = sum(xx[1] == 1)
      n_central_x2_0 = sum(xx[1] == 0)
      n_alter_x2_1 = sum(xx[-1] == 1)
      n_alter_x2_0 = sum(xx[-1] == 0) #these n's must sum to five
      
      p11 = (x2_prev + kappa)
      p01 = 1-p11
      p10 = x2_prev + (kappa / (length(aa) -1))
      p00 = 1 - p10
      
      p_x2_x1 = (p11^n_central_x2_1) * (p01^n_central_x2_0) * (p10^n_alter_x2_1) * (p00^n_alter_x2_0) 
      
      for (g in 1:ncol(gamma_numer)){
        #g = 1
        curr_gamma = gamma_numer[,g]
        
        lin_pred <- cbind(1, as.matrix(cbind(x1, xx))) %*% curr_gamma
        re_alpha <- FromAlphaToRE(alpha = p_trt, lin_pred = lin_pred,
                                  alpha_re_bound = 15)
        
        p_a_x = CalcNumerator(Ai_j = aa, #problem: these shoudl be probabilities, but some are greater than 1???
                              Xi_j = cbind(x1, xx),
                              gamma_numer = curr_gamma,
                              alpha = p_trt, 
                              re_alpha = re_alpha, 
                              include_alpha = FALSE)$prob
        res_mat[a,x,g] = p_a_x * p_x2_x1 #this is P(Aij | X2)
        
        #need outcome model to weight by this probability?
        
      }
    }
  }
  
  p_a_x1 = apply(res_mat, c(1,3), sum) #P(A | X1) #some of these are greater than 1
  
  wtd_y_components = p_a_x1 * as.vector(Y_mat[,1])
  wtd_y = apply(wtd_y_components, 2, sum)
  
  wtd_y0_components = p_a_x1 * as.vector(Y_mat[,2])
  wtd_y0 = apply(wtd_y0_components, 2, sum)
  
  wtd_y1_components = p_a_x1 * as.vector(Y_mat[,3])
  wtd_y1 = apply(wtd_y1_components, 2, sum)
  
  
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
#results = get_true_diffusion(gamma_numer = gamma_numer, x2_prev = 0.2, kappa = 0.25, trt_vecs = trt_vecs)

