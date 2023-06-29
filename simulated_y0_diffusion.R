#have Y_mat from the diffusion truth
#trt_vecs from the diffusion truth

##################################################
# Calculate ybar0 by simulation for each possible treatment vector
##################################################
trt_vecs = t(expand.grid(0:1, 0:1, 0:1, 0:1, 0:1))
out = rep(NA, 5)
p = 0.5 #p diffusion
out_tab = array(NA, c(500,5,32))

for (a in 1:ncol(trt_vecs)){
  #a = 5
  aa = as.vector(trt_vecs[,a])
  out = ifelse(aa==1, 1, NA)
  
  for(i in 1:500){
    if(aa[1] == 1){ #if center treated
      out[aa==0] = rbinom(length(out[aa==0]), 1, p)
    }
    
    if(aa[1] == 0){ #if center not treated
      if(sum(aa) == 0){out = 0}else{
        out[1] = max(rbinom(sum(aa[-1]), 1, p))
        out[2:5][aa[2:5]==0] = 0
      }
    }
    
    out = ifelse(aa==1, 1, out)
    
    out_tab[i,,a] = out
  }
}

avg = apply(out_tab, c(2,3), mean)

y0tab = c()
for (a in 1:ncol(trt_vecs)){
  aa = as.vector(trt_vecs[,a])
  yy = as.vector(avg[,a])
  y0tab = append(y0tab, mean((yy[aa == 0])))
}


#############################################
# SIMULATE 10,000 CLUSTERS WITH GAMMA = 0
#############################################

nclus = 10000
neigh = sort(rep(1:nclus, 5))
Aij = rbinom(50000, 1, 0.25)

clist = c()
for (cluster in 1:nclus){
  A = Aij[neigh == cluster]
  idx = which(sapply(1:ncol(trt_vecs), function(i) all(A==trt_vecs[,i]))==T)
  clist = append(clist, idx)
}
clist

###########################################
# WHAT IS THE TRUE YBAR0 OF GAMMA=0 SCENARIO
###########################################
mean(Y_mat[clist[clist!=32],2])
mean(y0tab[clist[clist!=32]])
