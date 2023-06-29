X1 <- c(1, 0, 0, 0, 0)
# probs <- c(0.33, 0.23, 0.23, 0.23, 0.23)
probs <- rep(0.25, 5)
mean(probs)
p_diff <- 0.5

B <- 5000
Ya <- array(NA, dim = c(4,B)) #Y0 of an alter
Yc <- rep(NA, B) #Y0 of a center unit

for (bb in 1 : B) {
  Ya_trt <- rbinom(4, 1, prob = probs[2:5]) #Draw treatment for the 4 alters (fix center trt=0)
  k <- sum(Ya_trt) #Number of treated alters
  if (k == 0) {
    Yc[bb] <- 0 #If no treated alters, center outcome is 0, no diffusion
  } else {
    Yc[bb] <- any(rbinom(k, 1, prob = p_diff) == 1) #If any alters are treated, outcome of center is bern(p_diff)
  }
  
  Yc_trt <- rbinom(1, 1, prob = probs[1]) #Draw treatment for the center (fix alter trt = 0)
  if (Yc_trt == 0) {
    Ya[1:4, bb] <- 0
  } else {
    Ya[1,bb] <- rbinom(1, 1, prob = p_diff)
    Ya[2,bb] <- rbinom(1, 1, prob = p_diff)
    Ya[3,bb] <- rbinom(1, 1, prob = p_diff)
    Ya[4,bb] <- rbinom(1, 1, prob = p_diff)
  }
}

clust_avg = rep(0, B)
for (bb in 1:B){
  clust_avg[bb] = (sum(Ya[,bb]) + Yc[bb]) / 5
}
mean(clust_avg)
mean(Yc)
mean(Ya)

Yave <- (4 / 5) * Ya + (1 / 5) * Yc
mean(Yave)

(1 - probs[1]) * 0 + probs[1] * p_diff
mean(Ya)

#analytical calculations

Ycent <- rep(0, 1000) #central unit
#colnames(Ycent) <- paste0('kk=', 0:4)
probs_k <- dbinom(1:4, 4, prob = probs[2])
for (kk in 1 : 4) {
  for (bb in 1 : 1000) {
    num_trt <- rbinom(1, 4, prob = probs_k)
    if (num_trt > 0) {
      Ycent[bb] <- any(rbinom(num_trt, 1, prob = p_diff) == 1)
    }
  }
}
mean(Ycent)
mean(Yc)



#######################################
#scenario: 
#one covar: center / alter indicator
#p_diffusion = 0.5
#p_trt / alpha = 0.25 for everyone


#Y0 of center: 
#p any diffuses to the center: 

#if all four are treated
1 - (1-diffusion_p)^4 #P(any diffuses to center is 0.9)

1 - (1-diffusion_p)^1 #with expectation of 1 treated alter

#p diff occurs
4*(.25/2) - 6*(.25/2)^2 - 3*(.25/2)^3 - (.25/2)

#p no diff occurs from unit = p unit not treated + p unit treated * diffusion doesn't occur
(1-p_trt) + (.25/2)


(.5*.25) + .75
