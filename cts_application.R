library(gamlss)
library(stringr)
library(tidyverse)

ngam = 11
gl = seq(from = -.2, to = .2, length.out = 10)
gamma_list = rbind(rep(0, 10),
                   c(gl),
                   rep(0, 10),
                   rep(0, 10),
                   rep(0, 10),
                   rep(0, 10))
gamma_cts = cbind(gamma_list, rep(0, 6))

load("/gpfs/gibbs/project/forastiere/sgd37/vax/vaccine_data_for_aim3.RSave")
df$neigh = as.numeric(factor(str_sub(df$census_geocode, 3, 5)))

