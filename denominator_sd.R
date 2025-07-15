# want options for probabilistic assignment and for n choose k assignment

#inputs: 
# A: trt vector for cluster
# X: covar matrix for cluster
# phi hat: the probability of treatment as observed in the data 
#        (num in treated in cluster / num in cluster)
# alpha: make up
# new: completely_random: is the data generated from a completely random experiment?
#                         default to no and use bernoulli. if yes, use 1 / binom coef

# Note: since the propensity is the same for all in cluster (bc indep of covars)
#       there should only be one denom value i think? or it makes sense to be the same 
#       only for given trt status

#test_trt = rbinom(20,1, .4)

calc_denominator_sd2 <- function(A, fix_phi = T, pop_p_trt, alpha) { #d i want to include dta as an input?
  
  #testing
  #A = df$trt[df$cluster == 2]
  #alpha = alpha
  #completely_random = FALSE
  #propensity_score = FALSE
  #dta = df
  #neigh_ind <- sapply(1 : max(dta$neigh), function(x) which(dta$neigh == x))
  #nn = 1
  #cov_cols = 33 
  
    n = length(A)
    k = sum((A==1))
    
    if(fix_phi == T){
      phi_hat = pop_p_trt#k/n #this is not phi_hat - p(trt) in cluster
    }else{
      phi_hat = k/n
    }
    
    
    p = ((phi_hat/alpha)^k) * ((1-phi_hat)/(1-alpha))^(n-k)
    #p = 1 / nchoosek(n,k) 
  
  
  return(p)
}

