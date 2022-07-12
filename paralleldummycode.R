# libraries ----
library(parallel)
library(pbmcapply)
library(ggplot2)
library(microbenchmark)

sammi_func <- function(clust_size){
  #for now clust size the same for every sim
  
  # a bunch of stuff happens that takes a long time
  
  #return a list of a bunch of arrays
  test1 = array(data = c(1:20), dim = c(2,10))
  test2 = array(data = c(20:40), dim = c(4,5))
  test3 = array(data = c(100:200), dim = c(10,10))
  
  
  return(list(array1 = test1,
              array2 = test2,
              array3 = test3,
              clust_size = clust_size))
}


# parallelize ----
n_sim <- 10
(n_cores <- detectCores())


#this works
output <- pbmclapply(1, mc.cores = n_cores, function(i) {
  res <- replicate(n_sim, sammi_func(i))
  res#c(mean(res), sd(res))
})


#this does not work!
# '15' is the argument that get_yhats_boot requires. it is called 'clust_size' in the function, like in the example
# I want to run identical function n_sim times with same params
output <- pbmclapply('15',mc.cores = n_cores, function(i) {
  res <- replicate(n_sim, get_yhats_boot(i))
  res#c(mean(res), sd(res))
}) #started at 12:48
Sys.time()
output[[1]][,1]

