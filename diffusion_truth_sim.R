library(here)
#new diffusion truth simulation
gl = seq(from = -.5, to = .5, length.out = 33)
ngl = rep(0, 33)
gamma_list = rbind(rep(0, 33),
                   c(ngl, gl),
                   c(gl, ngl))
gamma_list = cbind(gamma_list, c(0,0,0))
gamma_numer = gamma_list


#X1 = c(1,0,0,0,0)
n_sim = 2000
concors = c(0, 0.65)
diffusion_ps = c(0.2, 0.5, 0.8)
x2_prev = 0.2
n_sim_x2 = 200
clust_size = 5
alpha_re_bound = 20
alpha = 0.25


A_center = c(0, NA, NA, NA, NA)
A_alter = c(NA, 0)
Y0_truth = array(NA, dim = c(length(concors), length(diffusion_ps), n_sim_x2, n_sim, ncol(gamma_numer)))
X1 = c(1,0,0,0,0)
c = 0

for(concor in concors){
  c = c+1
  print(c)
  d = 0
for(diffusion_p in diffusion_ps){
  d = d+1
  print(d)
for(sim in 1:n_sim_x2){
  if(sim %/% 100){print(sim)}
  if(concor == 0){X2 = rbinom(clust_size,1, x2_prev)}
  if(concor == 0.65){
    kappa = 0.25
    X2 = rbinom(clust_size, 1, c(x2_prev + kappa, rep(x2_prev - kappa / (clust_size-1), clust_size-1)))
  }
  X = cbind(X1, X2)
  
  for (g in 1:ncol(gamma_numer)){
      #g = 60
    gamma_use = gamma_numer[,g]
    lin_pred <- cbind(1, as.matrix(X)) %*% gamma_use
      
    re_alpha <- FromAlphaToRE(alpha = alpha, lin_pred = lin_pred,
                                alpha_re_bound = alpha_re_bound)
      
    probs = expit(lin_pred + re_alpha)
    #mean(probs)
    
    for (cl in 1:n_sim){
      #cl =1
      ##### Yc0 #########
      new_Aj = rbinom(clust_size-1,1,probs[2:clust_size]) 
      #A_center[2:clust_size] = new_Aj
      k = sum(new_Aj) #number of alters treated
      diff_occurs = ifelse(k>0, max(rbinom(k, 1, diffusion_p)), 0)
      Yc0 = ifelse(diff_occurs == 1, 1, 0)
      
      ##### Ya0 #######
      A_alter[1] = rbinom(1,1,probs[1])
      diff_occurs = ifelse(A_alter[1] == 1, rbinom(1,1,diffusion_p), 0)
      Ya0 = ifelse(diff_occurs == 1, 1, 0)
      
      Ycl0 = (1/5)*Yc0 + (4/5)*Ya0 
      
      Y0_truth[c, d, sim, cl, g] <- Ycl0
      
  
    }
  }
}
}
}
#save(Y0_truth, file = here('Y0_truth_sim_concordance0.Rsave'))
#load(here('Y0_truth_sim_concordance0.Rsave'))
hist(apply(Y0_truth[1, 1, , , 60], 1, mean))
length(apply(Y0_truth[1, 1, , , 60], 2, mean))
str(Y0_truth)

Y0_truth_bar = apply(Y0_truth, c(1,2,5), mean)
dimnames(Y0_truth_bar) = list(concordance = c(0, 0.65), pdiff = c(0.2,0.5, 0.8), 
                              g = c(paste0(gamma_numer[3,1:33], '_X2'), paste0(gamma_numer[3,1:33], '_X1'), 
                                    NA))
Y0_truth_bar = as.data.frame.table(Y0_truth_bar) %>%
  extract(g, into = c("g", "gamma_ind"), "(.*)_([^_]+)$")
names(Y0_truth_bar)[5] <- 'true_y0'

Y0_truth_bar$true_y1 = 1
Y0_truth_bar$true_de = Y0_truth_bar$true_y1 - Y0_truth_bar$true_y0

Y0_truth_bar = Y0_truth_bar %>%
  merge(Y0_truth_bar %>% filter(is.na(gamma_ind)) %>% select(concordance, pdiff, ref_y0 = true_y0), by = c('concordance', 'pdiff'))

Y0_truth_bar$true_ie0 = Y0_truth_bar$true_y0 - Y0_truth_bar$ref_y0
Y0_truth_bar$true_ie1 = 0

save(Y0_truth_bar, file = here('Y0_truth_sim_df.Rsave'))


save.image('simtruthwkspc.Rsave')

load(here('simtruthwkspc.Rsave'))
load('/gpfs/ysm/project/forastiere/sgd37/cai/wkspcfeb28.Rsave')
