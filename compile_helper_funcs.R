#calc means over 600 simulations for each scenario
get_means = function(which_iter, e.names, beta4 = ugly_parm[[1]][3], diffusion = F){
  #which_iter = 1
  iters = alliters[str_detect(alliters, parmlist[which_iter])]
  
  for(i in e.names){
    assign(i, array(NA, dim = c(ngam, length(iters)), dimnames(list('gammas', 'iters'))))
  }
  
  #EXTRACT ALL THE ESTIMATES FROM THE ITERATIONS
  for (i in 1:length(iters)){
    #i = 1
    wk = iters[i]
    #wk = "scenario255_1_1_0.65.RSave"
    load(wk)
    
    #potential outcomes
    y0_ht[,i] = test$unadj[1,]
    y1_ht[,i] = test$unadj[2,]
    
    y0_haj[,i] = test$haj[1,]
    y1_haj[,i] = test$haj[2,]
    
    #direct effect and error
    # haj and ht
    direct_ht[,i] = test$ht_direct["est",]
    direct_analyticalvar_ht[,i] = test$ht_direct["var",]
    direct_bootvar_ht[,i] = test$ht_direct["boot_var",]
    
    direct_haj[,i] = test$direct["est",]
    direct_analyticalvar_haj[,i] = test$direct["var",]
    direct_bootvar_haj[,i] = test$direct["boot_var",]
    
    #indirect effect and error
    # haj and ht
    indirect0_haj[,i] = test$indirect0['est',]
    indirect0_analyticalvar_haj[,i] = test$indirect0['var',]
    indirect0_bootvar_haj[,i] = test$indirect0['boot_var',]
    
    indirect1_haj[,i] = test$indirect1['est',]
    indirect1_analyticalvar_haj[,i] = test$indirect1['var',]
    indirect1_bootvar_haj[,i] = test$indirect1['boot_var',]
    
    indirect0_ht[,i] = test$ht_indirect0['est',]
    indirect0_analyticalvar_ht[,i] = test$ht_indirect0['var',]
    indirect0_bootvar_ht[,i] = test$ht_indirect0['boot_var',]
    
    indirect1_ht[,i] = test$ht_indirect1['est',]
    indirect1_analyticalvar_ht[,i] = test$ht_indirect1['var',]
    indirect1_bootvar_ht[,i] = test$ht_indirect1['boot_var',]
    
    #overall effect and error
    #haj and ht
    oe_haj[,i] = test$oe['est',]
    oe_analyticalvar_haj[,i] = test$oe['var',]
    oe_bootvar_haj[,i] = test$oe['boot_var',]
    
    oe_ht[,i] = test$ht_oe['est',]
    oe_analyticalvar_ht[,i] = test$ht_oe['var',]
    oe_bootvar_ht[,i] = test$ht_oe['boot_var',]
    
    if(diffusion == F){
      true_ie0[,i] = test$het_ie_truth$ie[1,,ngam]
      true_ie1[,i] = test$het_ie_truth$ie[2,,ngam]
      true_de[,i] = 3
      true_oe[,i] = test$het_ie_truth$oe[,ngam]
      
    }else{
      true_ie0[,i] = test$truth$ie[1,,ngam]
      true_ie1[,i] = test$truth$ie[2,,ngam]
      true_de[,i] = test$truth$de
      true_oe[,i] = test$truth$oe[,ngam]
      true_y0[,i] = test$truth$y0[,ngam]
      true_y1[,i] = test$truth$y1[,ngam]
      
    }
    
    if(beta4 == 0){
      true_oe[,i] = 0
    }
    
    oe_mcvar_ht[,i] = test$ht_oe_mcvar
    
    #coverage
    oe_anacoverage_ht[,i] = (true_oe[,i] >= test$ht_oe['LB',]) & (true_oe[,i] <= test$ht_oe['UB',])
    oe_anacoverage_haj[,i] = (true_oe[,i] >= test$oe['LB',]) &  (true_oe[,i] <= test$oe['UB',])
    
    oe_bootcoverage_ht[,i] = (true_oe[,i] >= test$ht_oe['boot_var_LB',]) &  (true_oe[,i] <= test$ht_oe['boot_var_UB',]) 
    oe_bootcoverage_haj[,i] = (true_oe[,i] >= test$oe['boot_var_LB',]) & (true_oe[,i] <= test$oe['boot_var_UB',])
    
    #ie0 coverage
    ie0_anacoverage_ht[,i] = (true_ie0[,i] >= test$ht_indirect0['LB',]) & (true_ie0[,i] <= test$ht_indirect0['UB',])
    ie0_anacoverage_haj[,i] = (true_ie0[,i] >= test$indirect0['LB',]) &  (true_ie0[,i] <= test$indirect0['UB',])
    
    ie0_bootcoverage_ht[,i] = (true_ie0[,i] >= test$ht_indirect0['boot_var_LB',]) &  (true_ie0[,i] <= test$ht_indirect0['boot_var_UB',]) 
    ie0_bootcoverage_haj[,i] = (true_ie0[,i] >= test$indirect0['boot_var_LB',]) & (true_ie0[,i] <= test$indirect0['boot_var_UB',])
    
    
    #de coverage
    de_anacoverage_ht[,i] = (3 >= test$ht_direct['low_int',]) & (3 <= test$ht_direct['high_int',])
    de_anacoverage_haj[,i] = (3 >= test$direct['low_int',]) & (3 <= test$direct['high_int',])
    
    de_bootcoverage_ht[,i] = (3 >= test$ht_direct['boot_var_LB',]) & (3 <= test$ht_direct['boot_var_UB',])
    de_bootcoverage_haj[,i] = (3 >= test$direct['boot_var_LB',]) & (3 <= test$direct['boot_var_UB',])
    
    #ie coverage
    #ie_anacoverage_ht[,i] = (true_ie0[,i] > test$ht_oe['boot_var_LB',]) & (true_oe[,i] < test$ht_oe['UB',])
    #ie_anacoverage_haj[,i] = (true_ie0[,i] > test$oe['boot_var_LB',]) &  (true_oe[,i] < test$oe['UB',])
    
    #ie_bootcoverage_ht[,i] = (true_oe[,i] > test$ht_oe['boot_var_LB',]) &  (true_oe[,i] < test$ht_oe['boot_var_UB',]) 
    #ie_bootcoverage_haj[,i] = (true_oe[,i] > test$oe['boot_var_LB',]) & (true_oe[,i] < test$oe['boot_var_UB',])
  }
  
  #now i have the vector for each sim. next - average over the sims 43
  estimates = list(direct_ht, direct_analyticalvar_ht, direct_bootvar_ht, #HERE!!!
                   direct_haj, direct_analyticalvar_haj, direct_bootvar_haj,
                   indirect0_ht, indirect0_analyticalvar_ht, indirect0_bootvar_ht,
                   indirect0_haj, indirect0_analyticalvar_haj, indirect0_bootvar_haj,
                   indirect1_ht, indirect1_analyticalvar_ht, indirect1_bootvar_ht,
                   indirect1_haj, indirect1_analyticalvar_haj, indirect1_bootvar_haj,
                   oe_haj, oe_analyticalvar_haj, oe_bootvar_haj,
                   oe_ht, oe_analyticalvar_ht, oe_bootvar_ht,
                   true_ie0, true_ie1, true_oe, true_de, true_y0, true_y1,
                   y0_ht, y1_ht, y0_haj, y1_haj,
                   oe_anacoverage_ht, oe_anacoverage_haj, oe_bootcoverage_ht, oe_bootcoverage_haj,
                   de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj,
                   ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj,
                   oe_mcvar_ht
  )
  
  oe_mcvar_ht = apply(oe_ht, 1, var)
  
  #take the mean of the 600 simulations
  for(i in 1:length(e.names)){
    est = estimates[[i]]
    cal = apply(est, 1, mean, na.rm = T)
    if(e.names[i] != 'oe_mcvar_ht'){assign(e.names[i], cal)}
  }
  
  return(list(direct_ht, direct_analyticalvar_ht, direct_bootvar_ht, #HERE!!!
              direct_haj, direct_analyticalvar_haj, direct_bootvar_haj,
              indirect0_ht, indirect0_analyticalvar_ht, indirect0_bootvar_ht,
              indirect0_haj, indirect0_analyticalvar_haj, indirect0_bootvar_haj,
              indirect1_ht, indirect1_analyticalvar_ht, indirect1_bootvar_ht,
              indirect1_haj, indirect1_analyticalvar_haj, indirect1_bootvar_haj,
              oe_haj, oe_analyticalvar_haj, oe_bootvar_haj,
              oe_ht, oe_analyticalvar_ht, oe_bootvar_ht,
              true_ie0, true_ie1, true_oe, true_de, true_y0, true_y1,
              y0_ht, y1_ht, y0_haj, y1_haj,
              oe_anacoverage_ht, oe_anacoverage_haj, oe_bootcoverage_ht, oe_bootcoverage_haj,
              de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj,
              ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj,
              oe_mcvar_ht))
}



#pick wheterh to use boot or analytical CI's
pick_var <- function(v = c('a', 'b'), bigbiastab){
  if (v == 'b'){
    bigbiastab = bigbiastab %>% 
      mutate(de_ht_ana_lb = de_ht - 1.96*sqrt(direct_bootvar_ht),
             de_ht_ana_ub = de_ht + 1.96*sqrt(direct_bootvar_ht),
             ie0_ht_ana_lb = ie0_ht - 1.96*sqrt(indirect0_bootvar_ht),
             ie0_ht_ana_ub = ie0_ht + 1.96*sqrt(indirect0_bootvar_ht),
             ie1_ht_ana_lb = ie1_ht - 1.96*sqrt(indirect1_bootvar_ht),
             ie1_ht_ana_ub = ie1_ht + 1.96*sqrt(indirect1_bootvar_ht),
             oe_ht_ana_lb = oe_ht - 1.96*sqrt(oe_bootvar_ht),
             oe_ht_ana_ub = oe_ht + 1.96*sqrt(oe_bootvar_ht),
             de_haj_ana_lb = de_haj - 1.96*sqrt(direct_bootvar_haj),
             de_haj_ana_ub = de_haj + 1.96*sqrt(direct_bootvar_haj),
             ie0_haj_ana_lb = ie0_haj - 1.96*sqrt(indirect0_bootvar_haj),
             ie0_haj_ana_ub = ie0_haj + 1.96*sqrt(indirect0_bootvar_haj),
             ie1_haj_ana_lb = ie1_haj - 1.96*sqrt(indirect1_bootvar_haj),
             ie1_haj_ana_ub = ie1_haj + 1.96*sqrt(indirect1_bootvar_haj),
             oe_haj_ana_lb = oe_haj - 1.96*sqrt(oe_bootvar_haj),
             oe_haj_ana_ub = oe_haj + 1.96*sqrt(oe_bootvar_haj)) %>%
      filter(gamma_ind != 'X3')
  } else{
    bigbiastab = bigbiastab %>% 
      mutate(de_ht_ana_lb = de_ht - 1.96*sqrt(direct_analyticalvar_ht),
             de_ht_ana_ub = de_ht + 1.96*sqrt(direct_analyticalvar_ht),
             ie0_ht_ana_lb = ie0_ht - 1.96*sqrt(indirect0_analyticalvar_ht),
             ie0_ht_ana_ub = ie0_ht + 1.96*sqrt(indirect0_analyticalvar_ht),
             ie1_ht_ana_lb = ie1_ht - 1.96*sqrt(indirect1_analyticalvar_ht),
             ie1_ht_ana_ub = ie1_ht + 1.96*sqrt(indirect1_analyticalvar_ht),
             oe_ht_ana_lb = oe_ht - 1.96*sqrt(oe_analyticalvar_ht),
             oe_ht_ana_ub = oe_ht + 1.96*sqrt(oe_analyticalvar_ht),
             de_haj_ana_lb = oe_haj - 1.96*sqrt(direct_analyticalvar_haj),
             de_haj_ana_ub = oe_haj + 1.96*sqrt(direct_analyticalvar_haj),
             ie0_haj_ana_lb = ie0_haj - 1.96*sqrt(indirect0_analyticalvar_haj),
             ie0_haj_ana_ub = ie0_haj + 1.96*sqrt(indirect0_analyticalvar_haj),
             ie1_haj_ana_lb = ie1_haj - 1.96*sqrt(indirect1_analyticalvar_haj),
             ie1_haj_ana_ub = ie1_haj + 1.96*sqrt(indirect1_analyticalvar_haj),
             oe_haj_ana_lb = oe_haj - 1.96*sqrt(oe_analyticalvar_haj),
             oe_haj_ana_ub = oe_haj + 1.96*sqrt(oe_analyticalvar_haj)) %>%
      filter(gamma_ind != 'X3')
  }
  
  return(bigbiastab)
}

