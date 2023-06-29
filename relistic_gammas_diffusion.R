#CHECKING PLAUSIBLE GAMMAS FOR THE DIFFUSION SCENARIO

#GENERATE GIANT DATA SET
n_clus = 1000
clust_size = 5
neigh = sort(rep(1:n_clus, clust_size))
X1ij = rep(c(1,0,0,0,0), n_clus)
Aij = rbinom(n_clus*clust_size, 1, p_trt)
Yij = rep(NA, n_clus*clust_size)

sim = data.frame(neigh, X1ij, Aij, Yij)

#generate outcome
sim$is_center_trted <- rep(sim %>% group_by(neigh) %>% slice(1) %>% pull(Aij), each = 5)
sim$Tij <- rep(sim %>% group_by(neigh) %>% summarise(Tij = sum(Aij)) %>% pull(Tij), each = 5)

sim = sim %>%
  group_by(neigh) %>% 
  mutate(n_untrt_alter = sum(Aij==0 & X1ij==0 & is_center_trted == 1)) %>%
  mutate(Yij = case_when(Aij == 1 ~ 1, #treated node
                         #untrt alter where center IS NOT treated
                         Aij == 0 & X1ij == 0 & is_center_trted == 0 ~ 0,
                         #untreated center - #max of a coin flip for each treated alter in the cluster
                         Aij == 0 & X1ij == 1 ~ as.numeric(max(rbinom(Tij, 1, diffusion_p))), 
                         #untrt alter where center IS treated
                         Aij == 0 & X1ij == 0 & is_center_trted == 1 ~ as.numeric(rbinom(n_untrt_alter, 1, diffusion_p)))) #change the first 1 to number of rows

#generate x2 (indep)
sim$X2ij = rbinom(clust_size * n_clus, 1, x2_prev) #kappa = 0

dta <- sim %>% select(neigh, Aij, X1ij, X2ij, Yij, Tij)

c_list = c()
for(clus in unique(dta$neigh)){
  #clus = 1
  mini = dta %>% filter(neigh == clus)
  mod = glm(data = mini,
      Aij ~ X1ij,
      family = 'binomial')
  c_list = append(c_list, coef(mod)[2])
}
summary(c_list)
sd(c_list)
hist(c_list)
