make_simulated_dataset <- function(n_clus, n_i, n_time, n_gamma, concor,
                                   beta_0, beta_1, beta_2, beta_3, beta_4, beta_5, gamma_cts){
  
  sim_df = data.frame(neigh = c(), 
                      A = c(), 
                      X1 = c(), 
                      X2 = c(),
                      U = c(), V =c(), W =c(), 
                      time = c())
  
  for (neigh in 1:n_clus){
    #set treatment
    A_t1 = rbeta(n_i, .4, 1)
    A_t2 = expit(A_t1 + rnorm(n_i, 0,1))
    A_t0 = expit(A_t1 - rnorm(n_i, 0, 1))
    
    #set covars - not time varying
    X1 = rbinom(n_i, 1, 0.5)
    X2 = rbinom(n_i, 1, 0.5)
    if(concor == .65){
      X2[X1 == 1 & X2 == 0] <- rbinom(length(X2[X1 == 1 & X2 == 0]), 1, .255)
      X2[X1 == 0 & X2 == 1] <- rbinom(length(X2[X1 == 0 & X2 == 1]), 1, .745)
    }
    
    #set interference terms
    U_t1 = c(); for(ii in 1:n_i){U_t1 = append(U_t1, mean(A_t1[-ii]))}
    V_t1 = c(); for(ii in 1:n_i){V_t1 = append(V_t1, mean(A_t1[-ii][X1[-ii] == 1]))}
    W_t1 = c(); for(ii in 1:n_i){W_t1 = append(W_t1, mean(A_t1[-ii][X2[-ii] == 1]))}
    
    U_t2 = c(); for(ii in 1:n_i){U_t2 = append(U_t2, mean(A_t2[-ii]))}
    V_t2 = c(); for(ii in 1:n_i){V_t2 = append(V_t2, mean(A_t2[-ii][X1[-ii] == 1]))}
    W_t2 = c(); for(ii in 1:n_i){W_t2 = append(W_t2, mean(A_t2[-ii][X2[-ii] == 1]))}
    
    sim_df = rbind(sim_df, 
                   data.frame(neigh = rep(neigh, n_i),
                              A = c(A_t1, A_t2),
                              Alag = c(A_t0, A_t1),
                              X1 = rep(X1, 2),
                              X2 = rep(X2, 2),
                              U = c(U_t1, U_t2),
                              V = c(V_t1, V_t2), 
                              W = c(W_t1, W_t2), 
                              time = c(rep(1, n_i), rep(2, n_i))))
    sim_df = sim_df %>% mutate(V = ifelse(is.na(V), 0, V),
                               W = ifelse(is.na(W), 0, W))
  }
  sim_df$id = rep(1:(neigh*n_i), 2)
  Y = beta_0 + beta_1*sim_df$A + beta_2*sim_df$X1 + beta_3*sim_df$U + beta_4*sim_df$V + beta_5*sim_df$W
  sim_df$Y = exp(rnorm(nrow(sim_df), Y, 1))
  return(sim_df)
}
