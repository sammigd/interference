#calc means over 600 simulations for each scenario
get_means = function(which_iter, e.names, beta4 = b4, beta5 = b5, diffusion = F, truth_use = NULL){
  iters = alliters[str_detect(alliters, parmlist[which_iter])]
  print(ngam)
  for(i in e.names){
    assign(i, array(NA, dim = c(ngam, length(iters)), dimnames(list('gammas', 'iters'))))
  }

  #EXTRACT ALL THE ESTIMATES FROM THE ITERATIONS
  for (i in 1:length(iters)){
    #i = 1
    wk = iters[i]
    #wk = "scenario255_1_1_0.65.RSave"
    load(wk)
    
    ###############################
    # POTENTIAL OUTCOMES ##########
    ###############################
    #HT
    #print(str(y0_ht))
    #print(str(test$unadj))
    y0_ht[,i] = test$unadj[1,]
    y1_ht[,i] = test$unadj[2,]
    
    #HAJ
    y0_haj[,i] = test$haj[1,]
    y1_haj[,i] = test$haj[2,]
    
    #HT Y0 VAR
    y0_ht_var[,i] = test$ht_yhat_var[1,1,]
    y0_ht_bootvar[,i] = test$ypop_bootvar[1,]
    y0_ht_mcvar[,i] = y0_ht[,i]
    
    #HT Y1 VAR
    y1_ht_var[,i] = test$ht_yhat_var[2,2,]
    y1_ht_bootvar[,i] = test$ypop_bootvar[2,]
    y1_ht_mcvar[,i] = y1_ht[,i]
    
    #HT Y COV
    y_ht_covar[,i] = test$ht_yhat_var[1,2,]
    
    
    
    ###############################
    # DIRECT EFFECT ###############
    ###############################    
    
    #HT
    direct_ht[,i] = test$ht_direct["est",]
    direct_analyticalvar_ht[,i] = test$ht_direct["var",]
    direct_bootvar_ht[,i] = test$ht_direct["boot_var",]
    
    #HAJ
    direct_haj[,i] = test$direct["est",]
    direct_analyticalvar_haj[,i] = test$direct["var",]
    direct_bootvar_haj[,i] = test$direct["boot_var",]
    
    ###############################
    # INDIRECT EFFECT #############
    ###############################
    
    # HT
    indirect0_ht[,i] = test$ht_indirect0['est',]
    indirect0_analyticalvar_ht[,i] = test$ht_indirect0['var',]
    indirect0_bootvar_ht[,i] = test$ht_indirect0['boot_var',]
    
    indirect1_ht[,i] = test$ht_indirect1['est',]
    indirect1_analyticalvar_ht[,i] = test$ht_indirect1['var',]
    indirect1_bootvar_ht[,i] = test$ht_indirect1['boot_var',]
    
    # HAJ
    indirect0_haj[,i] = test$indirect0['est',]
    indirect0_analyticalvar_haj[,i] = test$indirect0['var',]
    indirect0_bootvar_haj[,i] = test$indirect0['boot_var',]
    
    indirect1_haj[,i] = test$indirect1['est',]
    indirect1_analyticalvar_haj[,i] = test$indirect1['var',]
    indirect1_bootvar_haj[,i] = test$indirect1['boot_var',]
    
    ###############################
    # OVERALL EFFECT ##############
    ###############################

    #HT
    oe_ht[,i] = test$ht_oe['est',]
    oe_analyticalvar_ht[,i] = test$ht_oe['var',]
    oe_bootvar_ht[,i] = test$ht_oe['boot_var',]
    
    #HAJ
    oe_haj[,i] = test$oe['est',]
    oe_analyticalvar_haj[,i] = test$oe['var',]
    oe_bootvar_haj[,i] = test$oe['boot_var',]
    
    #SUM OF WEGITHS
    sum_wts0[,i] = test$wtlist[1,]
    sum_wts1[,i] = test$wtlist[2,]
    
    
    
    ###############################
    # TRUE PARAMETERS #############
    ###############################
    
    if(diffusion == F){
      true_de[,i] = apply(test$truth$true_y_ie[,,2] - test$truth$true_y_ie[,,1], 2, mean, na.rm = T)
      true_y0[,i] = apply(test$truth$true_y_ie[,,1], 2, mean, na.rm = T)
      true_y1[,i] = apply(test$truth$true_y_ie[,,2], 2, mean, na.rm =T)
      
      true_ie0[,i] = test$truth$ie[1,,ngam]
      true_ie1[,i] = test$truth$ie[2,,ngam]
      true_oe[,i] = test$truth$oe[,ngam]
      
      if((beta4 == 0 | beta4 == '0') & beta5 == 0){
        true_oe[,i] = 0
        true_ie1[,i] = 0
        true_ie0[,i] = 0
      }
    }else{
      #true_de[,i] = truth_use$true_de#test$truth$de
      true_y0[,i] = NA#truth_use$true_y0#test$truth$y0
      true_y1[,i] = NA#truth_use$true_y1#test$truth$y1
      true_de[,i] = NA#truth_use$true_y1 - truth_use$true_y0
      
      true_ie0[,i] = NA#truth_use$true_ie0
      true_ie1[,i] = NA#truth_use$true_ie1
      true_oe[,i] = test$truth$oe[,67] #NA#truth_use$true_oe #NA
      
    }

    
    
    ###############################
    # COVERAGE ####################
    ###############################
    oe_anacoverage_ht[,i] = (true_oe[,i] >= test$ht_oe['LB',]) & (true_oe[,i] <= test$ht_oe['UB',])
    oe_anacoverage_haj[,i] = (true_oe[,i] >= test$oe['LB',]) &  (true_oe[,i] <= test$oe['UB',])
    
    oe_bootcoverage_ht[,i] = (true_oe[,i] >= test$ht_oe['boot_var_LB',]) &  (true_oe[,i] <= test$ht_oe['boot_var_UB',]) 
    oe_bootcoverage_haj[,i] = (true_oe[,i] >= test$oe['boot_var_LB',]) & (true_oe[,i] <= test$oe['boot_var_UB',])
    
    #ie0 coverage
    ie0_anacoverage_ht[,i] = (true_ie0[,i] >= test$ht_indirect0['LB',]) & (true_ie0[,i] <= test$ht_indirect0['UB',])
    ie0_anacoverage_haj[,i] = (true_ie0[,i] >= test$indirect0['LB',]) &  (true_ie0[,i] <= test$indirect0['UB',])
    
    ie0_bootcoverage_ht[,i] = (true_ie0[,i] >= test$ht_indirect0['boot_var_LB',]) &  (true_ie0[,i] <= test$ht_indirect0['boot_var_UB',]) 
    ie0_bootcoverage_haj[,i] = (true_ie0[,i] >= test$indirect0['boot_var_LB',]) & (true_ie0[,i] <= test$indirect0['boot_var_UB',])
    
    #ie1 coverage
    ie1_anacoverage_ht[,i] = (true_ie1[,i] >= test$ht_indirect1['LB',]) & (true_ie1[,i] <= test$ht_indirect1['UB',])
    ie1_anacoverage_haj[,i] = (true_ie1[,i] >= test$indirect1['LB',]) &  (true_ie1[,i] <= test$indirect1['UB',])
    
    ie1_bootcoverage_ht[,i] = (true_ie1[,i] >= test$ht_indirect1['boot_var_LB',]) &  (true_ie1[,i] <= test$ht_indirect1['boot_var_UB',]) 
    ie1_bootcoverage_haj[,i] = (true_ie1[,i] >= test$indirect1['boot_var_LB',]) & (true_ie1[,i] <= test$indirect1['boot_var_UB',])
    
    #de coverage
    de_anacoverage_ht[,i] = (true_de[,i] >= test$ht_direct['low_int',]) & (true_de[,i] <= test$ht_direct['high_int',])
    de_anacoverage_haj[,i] = (true_de[,i] >= test$direct['low_int',]) & (true_de[,i] <= test$direct['high_int',])
    
    de_bootcoverage_ht[,i] = (true_de[,i] >= test$ht_direct['boot_var_LB',]) & (true_de[,i] <= test$ht_direct['boot_var_UB',])
    de_bootcoverage_haj[,i] = (true_de[,i] >= test$direct['boot_var_LB',]) & (true_de[,i] <= test$direct['boot_var_UB',])
    
    #y0 coverage
    y0_anacoverage_ht[,i] = (true_y0[,i] >= y0_ht[,i] - 1.96*sqrt(y0_ht_var[,i])) & 
      (true_y0[,i] <= y0_ht[,i] + 1.96*sqrt(y0_ht_var[,i]))
    #y1 coverage
    y1_anacoverage_ht[,i] = (true_y1[,i] >= y1_ht[,i] - 1.96*sqrt(y1_ht_var[,i])) & 
      (true_y1[,i] <= y1_ht[,i] + 1.96*sqrt(y1_ht_var[,i]))
    #y0 boot coverage
    y0_bootcoverage_ht[,i] = (true_y0[,i] >= y0_ht[,i] - 1.96*sqrt(y0_ht_bootvar[,i])) & 
      (true_y0[,i] <= y0_ht[,i] + 1.96*sqrt(y0_ht_bootvar[,i]))
    #y1 boot coverage
    y1_bootcoverage_ht[,i] = (true_y1[,i] >= y1_ht[,i] - 1.96*sqrt(y1_ht_bootvar[,i])) & 
      (true_y1[,i] <= y1_ht[,i] + 1.96*sqrt(y1_ht_bootvar[,i]))
    
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
                   y0_ht_var, y1_ht_var, y0_ht_bootvar, y1_ht_bootvar, y0_ht_mcvar, y1_ht_mcvar, 
                   oe_anacoverage_ht, oe_anacoverage_haj, oe_bootcoverage_ht, oe_bootcoverage_haj,
                   de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj,
                   ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj,
                   ie1_anacoverage_ht, ie1_anacoverage_haj, ie1_bootcoverage_ht, ie1_bootcoverage_haj,
                   y0_anacoverage_ht, y1_anacoverage_ht, y_ht_covar, y0_bootcoverage_ht, y1_bootcoverage_ht, sum_wts0, sum_wts1
  )
  
  
  #take the mean of the 600 simulations
  for(i in 1:length(e.names)){
    #print(str(estimates[[i]]))
    est = estimates[[i]]
    if(e.names[i] == 'y0_ht_mcvar' | e.names[i] == 'y1_ht_mcvar'){cal = apply(est, 1, var)}else{cal = apply(est, 1, mean)}
    assign(e.names[i], cal)
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
              y0_ht_var, y1_ht_var, y0_ht_bootvar, y1_ht_bootvar, y0_ht_mcvar, y1_ht_mcvar, 
              oe_anacoverage_ht, oe_anacoverage_haj, oe_bootcoverage_ht, oe_bootcoverage_haj,
              de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj,
              ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj,
              ie1_anacoverage_ht, ie1_anacoverage_haj, ie1_bootcoverage_ht, ie1_bootcoverage_haj,
              y0_anacoverage_ht, y1_anacoverage_ht, y_ht_covar, y0_bootcoverage_ht, y1_bootcoverage_ht, 
              sum_wts0, sum_wts1))
}



#pick wheterh to use boot or analytical CI's
pick_var <- function(v = c('a', 'b'), bigbiastab){
  if (v == 'b'){
    bigbiastab = bigbiastab %>% 
      mutate(de_ht_ana_lb = de_ht - 1.96*sqrt(de_bootvar_ht),
             de_ht_ana_ub = de_ht + 1.96*sqrt(de_bootvar_ht),
             ie0_ht_ana_lb = ie0_ht - 1.96*sqrt(indirect0_bootvar_ht),
             ie0_ht_ana_ub = ie0_ht + 1.96*sqrt(indirect0_bootvar_ht),
             ie1_ht_ana_lb = ie1_ht - 1.96*sqrt(indirect1_bootvar_ht),
             ie1_ht_ana_ub = ie1_ht + 1.96*sqrt(indirect1_bootvar_ht),
             oe_ht_ana_lb = oe_ht - 1.96*sqrt(oe_bootvar_ht),
             oe_ht_ana_ub = oe_ht + 1.96*sqrt(oe_bootvar_ht),
             de_haj_ana_lb = de_haj - 1.96*sqrt(de_bootvar_haj),
             de_haj_ana_ub = de_haj + 1.96*sqrt(de_bootvar_haj),
             ie0_haj_ana_lb = ie0_haj - 1.96*sqrt(indirect0_bootvar_haj),
             ie0_haj_ana_ub = ie0_haj + 1.96*sqrt(indirect0_bootvar_haj),
             ie1_haj_ana_lb = ie1_haj - 1.96*sqrt(indirect1_bootvar_haj),
             ie1_haj_ana_ub = ie1_haj + 1.96*sqrt(indirect1_bootvar_haj),
             oe_haj_ana_lb = oe_haj - 1.96*sqrt(oe_bootvar_haj),
             oe_haj_ana_ub = oe_haj + 1.96*sqrt(oe_bootvar_haj)) #%>%
      #filter(gamma_ind != 'X3')
  } else{
    bigbiastab = bigbiastab %>% 
      mutate(de_ht_ana_lb = de_ht - 1.96*sqrt(de_analyticalvar_ht),
             de_ht_ana_ub = de_ht + 1.96*sqrt(de_analyticalvar_ht),
             ie0_ht_ana_lb = ie0_ht - 1.96*sqrt(indirect0_analyticalvar_ht),
             ie0_ht_ana_ub = ie0_ht + 1.96*sqrt(indirect0_analyticalvar_ht),
             ie1_ht_ana_lb = ie1_ht - 1.96*sqrt(indirect1_analyticalvar_ht),
             ie1_ht_ana_ub = ie1_ht + 1.96*sqrt(indirect1_analyticalvar_ht),
             oe_ht_ana_lb = oe_ht - 1.96*sqrt(oe_analyticalvar_ht),
             oe_ht_ana_ub = oe_ht + 1.96*sqrt(oe_analyticalvar_ht),
             de_haj_ana_lb = de_haj - 1.96*sqrt(de_analyticalvar_haj),
             de_haj_ana_ub = de_haj + 1.96*sqrt(de_analyticalvar_haj),
             ie0_haj_ana_lb = ie0_haj - 1.96*sqrt(indirect0_analyticalvar_haj),
             ie0_haj_ana_ub = ie0_haj + 1.96*sqrt(indirect0_analyticalvar_haj),
             ie1_haj_ana_lb = ie1_haj - 1.96*sqrt(indirect1_analyticalvar_haj),
             ie1_haj_ana_ub = ie1_haj + 1.96*sqrt(indirect1_analyticalvar_haj),
             oe_haj_ana_lb = oe_haj - 1.96*sqrt(oe_analyticalvar_haj),
             oe_haj_ana_ub = oe_haj + 1.96*sqrt(oe_analyticalvar_haj)) #%>%
      #filter(gamma_ind != 'X3')
  }
  
  return(bigbiastab)
}


make_ms_figtab <- function(bigbiastab, effect, diffusion = F){
  #if(effect == 'OE'){bigbiastab = bigbiastab %>% pivot_longer(c(true_oe, oe_haj)); ylab = 'Overall Effect'}
  #if(effect == 'DE'){bigbiastab = bigbiastab %>% pivot_longer(c(true_de, de_haj)); ylab = 'Direct Effect'}
  #if(effect == 'IE0'){bigbiastab = bigbiastab %>% pivot_longer(c(true_ie0, ie0_haj)); ylab = 'Indirect Effect (0)'}
  #if(effect == 'IE1'){bigbiastab = bigbiastab %>% pivot_longer(c(true_ie1, ie1_haj)); ylab = 'Indirect Effect (1)'}
  
  lbx = bigbiastab %>% pull(names(bigbiastab)[str_detect(names(bigbiastab), paste0(str_to_lower(effect), '_haj_ana_lb'))])
  ubx = bigbiastab %>% pull(names(bigbiastab)[str_detect(names(bigbiastab), paste0(str_to_lower(effect), '_haj_ana_ub'))])
  
  fig2df = bigbiastab %>%
    mutate(lb = lbx, ub = ubx)
  
  a = 'No Interference\n (\U03B2' %p% subsc('3') %p% '=0; \U03B2' %p% subsc('4') %p% '=0)'
  b = 'Homogeneous \nInterference\n (\U03B2' %p% subsc('3') %p% '=1; \U03B2' %p% subsc('4') %p% '=0)'
  c = 'Moderate \nHeterogeneous \nInterference\n (\U03B2' %p% subsc('3') %p% '=0; \U03B2' %p% subsc('4') %p% '=1)'
  d = 'Strong \nHeterogeneous \nInterference\n (\U03B2' %p% subsc('3') %p% '=0; \U03B2' %p% subsc('4') %p% '=2)'

  e = 'No Interference\n (\U03B2' %p% subsc('5') %p% '=0)'
  f = 'Heterogeneous Interference\n (\U03B2' %p% subsc('5') %p% '=1)'
  
  if(!diffusion){
    fig2df = fig2df %>% 
      filter(b3 == 0  | 
             b3 == 1 & b4 == 0) %>%
      mutate(Interference = case_when(b3 == 0 & b4 == 0  ~ a,
                                      b3 == 1 & b4 == 0  ~ b,
                                      b3 == 0 &  b4 == 1 ~ c,
                                      b3 == 0 & b4 == 2  ~ d)) %>%  
      filter(!is.na(Interference)) %>%
      mutate(Interference = factor(Interference, levels = c(a,
                                                            b,
                                                            c,
                                                            d))) %>%
      mutate(InterferenceX2 = ifelse(b5 == 0, e, f),
             InterferenceX2 = factor(InterferenceX2, levels = c(e, 
                                                               f)))
      
  }
  if(diffusion){
    fig2df = fig2df %>%
      mutate(pdiff_long = paste('P(diffusion) = ', pdiff))
  }
  #print(table(fig2df$concordance))
  a = 'Cor(X' %p% supsc('(1)') %p%',X' %p% supsc('(2)') %p% ')=0'
  b = 'Cor(X' %p% supsc('(1)') %p%',X' %p% supsc('(2)') %p% ')!=0'
  fig2df = fig2df %>%
    mutate(Conc = case_when(concordance == '0' | concordance == 0 | concordance == '0.RSave' ~ a,
                            concordance == 0.65 | concordance == '0.8' | concordance == '0.65' | concordance == '0.65.RSave' ~ b)) %>%
    mutate(Conc = factor(Conc, levels = c(a,b))) %>%#,
           #name = ifelse(str_detect(name, 'haj'), paste0('Hajek ', effect), paste0('True ', effect)),
           #name = factor(name, levels = c( paste0('True ', effect),paste0('Hajek ', effect)))) %>%
    mutate(gamma_ind = case_when(gamma_ind == 'X1' ~ 'Target X' %p% supsc('(1)'),
                                 gamma_ind == 'X2' ~ 'Target X' %p% supsc('(2)')))
  return(list(fig2df = fig2df, effect = effect))
}

FromAlphaToRE <- function(alpha, lin_pred, alpha_re_bound = 10) {
  
  alpha_re_bound <- abs(alpha_re_bound)
  
  r <- optimise(f = AlphaToBi, lower = - 1 * alpha_re_bound,
                upper = alpha_re_bound, alpha = alpha, lin_pred = lin_pred)
  r <- r$minimum
  if (alpha_re_bound - abs(r) < 0.1) {
    warning(paste0('bi = ', r, ', alpha_re_bound = ', alpha_re_bound))
  }
  return(r)
}

AlphaToBi <- function(b, alpha, lin_pred) {
  exp_lin_pred <- exp(lin_pred)
  r <- abs(mean(exp_lin_pred / (exp_lin_pred + exp(- b))) - alpha)
  return(r)
}
