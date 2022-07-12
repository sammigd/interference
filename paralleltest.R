# libraries ----
library(parallel)
library(pbmcapply)
library(microbenchmark)

source('creatingfunction.R')

n_sim <- 10
(n_cores <- detectCores())


#this does not work on hpc
output <- pbmclapply('15', mc.cores = n_cores, function(i) {
  res <- replicate(n_sim, get_yhats_boot(i))
  res
}) 

#this works on hpc
output <- lapply('15', function(i) {
  res <- replicate(n_sim, get_yhats_boot(i))
  res
}) 

