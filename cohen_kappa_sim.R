X1ij = rbinom(1000,1, 0.5)

X2ij = rbinom(1000,1, 0.5)

if(concor == 0.65){ 
  X2ij[X1ij == 1 & X2ij == 0] <- rbinom(length(X2ij[X1ij == 1 & X2ij == 0]), 1, .255)
  X2ij[X1ij == 0 & X2ij == 1] <- rbinom(length(X2ij[X1ij == 0 & X2ij == 1]), 1, .745)
}

if(concor == 0.8){
  X2ij[X1ij == 1 & X2ij == 0] <- rbinom(length(X2ij[X1ij == 1 & X2ij == 0]), 1, .55)
  X2ij[X1ij == 0 & X2ij == 1] <- rbinom(length(X2ij[X1ij == 0 & X2ij == 1]), 1, .45)
}

#calculate cohen's kappa
res = mean(X1ij == X2ij)
pyes <- sum(X1ij) * sum(X2ij) / (1000*1000)
pno <- sum(1 - X1ij) * sum(1 - X1ij) / (1000*1000)

pe <- pyes + pno
coh <- (res - pe) / (1 - pe)
coh


#if i want cohens kappa to be approx 0.27 
#(correspeonding to 0.65 concordance for the non-diffusion scenario)

#fixed: prevalence, cohen's kapp

#vary: kappa?
prev_x2 = rep(NA, 1000)
kappa <- 0#0.25
prev <- 0.2
for (bb in 1: 1000) {
  x1 <- c(1, 0, 0, 0, 0)
  x2 <- rbinom(5, 1, c(prev + kappa, rep(prev - kappa / 4, 4)))
  res[bb] <- mean(x1 == x2)
  prev_x2[bb] <- mean(x2)
  
  pyes <- sum(x1) * sum(x2) / 25
  pno <- sum(1 - x1) * sum(1 - x2) / 25
  
  pe[bb] <- pyes + pno
  coh[bb] <- (res[bb] - pe[bb]) / (1 - pe[bb])
}
mean(res)
mean(prev_x2)
mean(pe)
mean(coh)
hist(coh)


calc_cohen_k<- function(x1, x2){
  res <- mean(x1 == x2)
  prev_x2 <- mean(x2)
  
  pyes <- sum(x1) * sum(x2) / length(x1)^2
  pno <- sum(1 - x1) * sum(1 - x2) / length(x1)^2
  
  pe <- pyes + pno
  coh <- (res - pe) / (1 - pe)
  return(coh)
}
