setwd("~/project/cai/diffusion_2023nov16")


alliters = list.files()
table(str_remove(alliters, word(alliters, sep = '_')))
parmlist = unique(str_remove(alliters, word(alliters, sep = '_')))#[-7]

gl = seq(from = -.2, to = .2, length.out = 33)
ngl = rep(0, 33)
gamma_list = rbind(rep(0, 33),
                   c(ngl, gl),
                   c(gl, ngl))
gamma_list = cbind(gamma_list, c(0,0,0))
gamma_numer = gamma_list
ngam = ncol(gamma_numer)

true_oe_df = data.frame(g = c(), gamma_ind = c(), concordance = c(), pdiff = c(),true_oe = c())

for (parms in 1:length(parmlist)){
  print(paste('parms=',parms))
  #make parms pretty
  ugly_parm = str_split(parmlist[parms], '_')
  pretty_parm = paste0('conc = ', ugly_parm[[1]][2], ', \np(diff) = ', str_sub(ugly_parm[[1]][3], start = 1, end =-7 ))
  concordance = ugly_parm[[1]][2]
  pdiff = str_sub(ugly_parm[[1]][3], start = 1, end =-7 )
  
  iters = alliters[str_detect(alliters, parmlist[parms])][1]
  
  load(iters)
  
  true_oe = test$truth$oe[,ngam]
  
  new_df = data.frame(g = c(rep(gamma_list[3,1:floor(ngam/2)], 2), 0),
                      gamma_ind = c(rep('X2', 33), rep('X1', 33), 'X1'),
                      concordance = concordance,
                      pdiff = pdiff,
                      true_oe = true_oe)
  true_oe_df = rbind(true_oe_df, new_df)
  
}
  