#library(purrr)
#library(abind)

setwd("~/project/cai/bigsims")

#final proj used the OG analytical variance for ht from GA code

iters = list.files()
iterlab = str_split(str_remove(iters, '.RSave'), '_')
iterlab = as.data.frame(abind(iterlab, along = 0))
iterlab$pgroup = paste0(iterlab$V2, iterlab$V3, iterlab$V4)
iterlab %>% group_by(pgroup) %>% summarise(n())

iters = iters[str_detect(iters, '_0_1_0.RSave')]


#CREATE EMPTY ARRAYS FOR EACH OF THE THINGS WE ARE ESTIMATING
e.names = c('direct_ht', 'direct_analyticalvar_ht', 'direct_bootvar_ht',
            'direct_haj', 'direct_analyticalvar_haj', 'direct_bootvar_haj',
            'indirect0_ht', 'indirect0_analyticalvar_ht', 'indirect0_bootvar_ht',
            'indirect0_haj', 'indirect0_analyticalvar_haj', 'indirect0_bootvar_haj',
            'indirect1_ht', 'indirect1_analyticalvar_ht', 'indirect1_bootvar_ht',
            'indirect1_haj', 'indirect1_analyticalvar_haj', 'indirect1_bootvar_haj',
            'oe_haj', 'oe_analyticalvar_haj', 'oe_bootvar_haj',
            'oe_ht', 'oe_analyticalvar_ht', 'oe_bootvar_ht'
            )

for(i in e.names){
  assign(i, array(NA, dim = c(126, length(iters)), dimnames(list('gammas', 'iters'))))
}

#EXTRACT ALL THE ESTIMATES FROM THE ITERATIONS
for (i in 1:length(iters)){
  wk = iters[i]
  #wk = "nointerference1.RSave"
  load(wk)
  
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
  
  #yhats
}

#compare averages
#apply(oe_ht, 1, mean)
#apply(sqrt(oe_bootvar_ht), 1, mean)
#apply(sqrt(oe_analyticalvar_ht), 1, mean)



#calculate coverage for iebootht and ieanalytical ht
a = array(NA, dim = c(126, length(iters)))
ie0.cov.boot.ht = a
ie0.cov.boot.haj = a

ie0.cov.ana.ht = a
ie0.cov.ana.haj = a

de.cov.ana.ht = a
de.cov.boot.ht = a
de.cov.ana.haj = a
de.cov.boot.haj = a


oe.cov.boot.haj = a
oe.cov.boot.ht = a
oe.cov.ana.haj = a
oe.cov.ana.ht = a


for (i in 1:length(iters)){
  for (gam in 1:126){
    ie0.cov.boot.ht[gam,i] = (0 >= indirect0_ht[gam,i] - 1.96*sqrt(indirect0_bootvar_ht[gam,i]) &
                             0 <= indirect0_ht[gam,i] + 1.96*sqrt(indirect0_bootvar_ht[gam,i]))
    
    ie0.cov.boot.haj[gam,i] = (0 >= indirect0_haj[gam,i] - 1.96*sqrt(indirect0_bootvar_haj[gam,i]) &
                              0 <= indirect0_haj[gam,i] + 1.96*sqrt(indirect0_bootvar_haj[gam,i]))
    ie0.cov.ana.ht[gam,i] = between(0, 
                                 indirect0_ht[gam,i] - 1.96*sqrt(indirect0_analyticalvar_ht[gam,i]),
                                 indirect0_ht[gam,i] + 1.96*sqrt(indirect0_analyticalvar_ht[gam,i]))
    ie0.cov.ana.haj[gam,i] = between(0, 
                                 indirect0_haj[gam,i] - 1.96*sqrt(indirect0_analyticalvar_haj[gam,i]),
                                 indirect0_haj[gam,i] + 1.96*sqrt(indirect0_analyticalvar_haj[gam,i]))
    #
    de.cov.ana.ht[gam,i] = between(3, 
                                 direct_ht[gam,i] - 1.96*sqrt(direct_analyticalvar_ht[gam,i]), 
                                 direct_ht[gam,i] + 1.96*sqrt(direct_analyticalvar_ht[gam,i]))
    
    de.cov.boot.ht[gam,i] = between(3, 
                                direct_ht[gam,i] - 1.96*sqrt(direct_bootvar_ht[gam,i]), 
                                direct_ht[gam,i] + 1.96*sqrt(direct_bootvar_ht[gam,i]))
    
    de.cov.ana.haj[gam,i] = between(3, 
                                   direct_haj[gam,i] - 1.96*sqrt(direct_analyticalvar_haj[gam,i]), 
                                   direct_haj[gam,i] + 1.96*sqrt(direct_analyticalvar_haj[gam,i]))
    de.cov.boot.haj[gam,i] = between(3, 
                                    direct_haj[gam,i] - 1.96*sqrt(direct_bootvar_haj[gam,i]), 
                                    direct_haj[gam,i] + 1.96*sqrt(direct_bootvar_haj[gam,i]))
    
    oe.cov.boot.haj[gam,i] = between(0, 
                             oe_haj[gam,i] - 1.96*sqrt(oe_bootvar_haj[gam,i]), 
                             oe_haj[gam,i] + 1.96*sqrt(oe_bootvar_haj[gam,i]))
    oe.cov.boot.ht[gam,i] = between(0, 
                             oe_ht[gam,i] - 1.96*sqrt(oe_bootvar_ht[gam,i]), 
                             oe_ht[gam,i] + 1.96*sqrt(oe_bootvar_ht[gam,i]))
    oe.cov.ana.haj[gam,i] = between(0, 
                                  oe_haj[gam,i] - 1.96*sqrt(oe_analyticalvar_haj[gam,i]), 
                                  oe_haj[gam,i] + 1.96*sqrt(oe_analyticalvar_haj[gam,i]))
    oe.cov.ana.ht[gam,i] = between(0, 
                                 oe_ht[gam,i] - 1.96*sqrt(oe_analyticalvar_ht[gam,i]), 
                                 oe_ht[gam,i] + 1.96*sqrt(oe_analyticalvar_ht[gam,i]))
  }
}

ie0.cov.boot.ht = apply(ie0.cov.boot.ht, 1, mean) #ok
ie0.cov.boot.haj = apply(ie0.cov.boot.haj, 1, mean) #ok

ie0.cov.ana.ht = apply(ie0.cov.ana.ht, 1, mean) #ok
ie0.cov.ana.haj = apply(ie0.cov.ana.haj, 1, mean) #ok 

de.cov.ana.ht = apply(de.cov.ana.ht,1,mean)
de.cov.boot.ht = apply(de.cov.boot.ht,1,mean)

de.cov.ana.haj = apply(de.cov.ana.haj,1,mean)
de.cov.boot.haj = apply(de.cov.boot.haj,1,mean)

oe.cov.boot.haj = apply(oe.cov.boot.haj, 1, mean)
oe.cov.ana.haj = apply(oe.cov.ana.haj, 1, mean)
oe.cov.boot.ht = apply(oe.cov.boot.ht, 1, mean)
oe.cov.ana.ht = apply(oe.cov.ana.ht, 1, mean)

covtab = data.frame(gamma = rep(gamma_list[2,], 6),
                    type = c(rep('HT', 3*126), rep('Hajek', 3*126)),
                    effect = rep(c(rep('Direct Effect', 126),
                                    rep('Indirect Effect', 126), 
                                    rep('Overall Effect', 126)),2),
                    est = c(de.cov.boot.ht, ie0.cov.boot.ht, oe.cov.boot.ht,
                            de.cov.boot.haj,ie0.cov.boot.haj,oe.cov.boot.haj
                            ))
ggplot(data = covtab, 
       aes(x = gamma, y = est, group = type)) + 
  geom_line() + facet_wrap(~type + effect) + 
  theme(text = element_text(size = 12)) + 
  labs(title = 'Bootstrap variance coverages')

covtab2 = data.frame(gamma = rep(gamma_list[2,], 6),
                    type = c(rep('HT', 3*126), rep('Hajek', 3*126)),
                    effect = rep(c(rep('Direct Effect', 126),
                                   rep('Indirect Effect', 126), 
                                   rep('Overall Effect', 126)),2),
                    est = c(de.cov.ana.ht, ie0.cov.ana.ht, oe.cov.ana.ht,
                            de.cov.ana.haj,ie0.cov.ana.haj,oe.cov.ana.haj
                    ))
ggplot(data = covtab2, 
       aes(x = gamma, y = est, group = type)) + 
  geom_point() + facet_wrap(~type + effect) + 
  theme(text = element_text(size = 12)) + 
  geom_hline(yintercept = .95, col = 'red') + 
  labs('analytical variance')

allcovtab = rbind(covtab %>% mutate(variance = 'Bootstrap'),
                  covtab2 %>% mutate(variance = 'Analytical'))

ggplot(data = allcovtab,aes(x = gamma, y = est, colour = variance)) + 
  geom_line() + 
  facet_wrap(~type + effect) + 
  theme(text = element_text(size = 12)) + 
  geom_hline(yintercept = .95, col = 'pink') + 
  labs('coverages')


apply(track, 1, mean)
#AVERAGE OVER EACH SIMULATION
estimates = list(direct_ht, direct_analyticalvar_ht, direct_bootvar_ht, #HERE!!!
                 direct_haj, direct_analyticalvar_haj, direct_bootvar_haj,
                 indirect0_ht, indirect0_analyticalvar_ht, indirect0_bootvar_ht,
                 indirect0_haj, indirect0_analyticalvar_haj, indirect0_bootvar_haj,
                 indirect1_ht, indirect1_analyticalvar_ht, indirect1_bootvar_ht,
                 indirect1_haj, indirect1_analyticalvar_haj, indirect1_bootvar_haj,
                 oe_haj, oe_analyticalvar_haj, oe_bootvar_haj,
                 oe_ht, oe_analyticalvar_ht, oe_bootvar_ht)

for(i in 1:length(e.names)){
  est = estimates[[i]]
  cal = apply(est, 1, mean)
  assign(e.names[i], cal)
}

    
#DIRECT VARIANCE COMPARISON
devartab = data.frame(gamma = gamma_list[2,],
                      hajest = direct_haj,
                      hajboot = direct_bootvar_haj,
                      hajanalytical = direct_analyticalvar_haj,
                      htest = direct_ht,
                      htboot = direct_bootvar_ht,
                      htanalytical = direct_analyticalvar_ht) %>% pivot_longer(hajest:htanalytical)
ggplot(devartab %>% filter(!(name %in% c('hajest', 'htest'))), aes(x = gamma, y = value)) + 
  facet_wrap(~name) + 
  geom_line() 
str(devartab)

devartab1 = data.frame(gamma = gamma_list[2,],
                      est = direct_haj,
                      boot = direct_bootvar_haj,
                      analytical = direct_analyticalvar_haj,
                      type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance')
devartab2 = data.frame(gamma = gamma_list[2,],
                      est = direct_ht,
                      boot = direct_bootvar_ht,
                      analytical = direct_analyticalvar_ht,
                      type = 'Horvitz Thompson') %>% pivot_longer(boot:analytical, 'variance')
devartab = rbind(devartab1, devartab2) %>%
  mutate(lb = est - 1.96*sqrt(value),
         ub = est + 1.96*sqrt(value),
         variance = ifelse(variance == 'boot', 'Bootstrapped Variance', 'Analytical Variance'))

ggplot(devartab, aes(gamma, est)) + 
  geom_line() + facet_wrap(~type+variance) + 
  geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), fill="blue", alpha=0.25) + 
  labs(y = 'direct effect',
       title = 'Direct Effect Estimators and 95% Confidence Intervals') + 
  theme(text = element_text(size = 14))

#INDIRECT VARIANCE COMPARISON
ievartab = data.frame(gamma = gamma_list[2,],
                      #ie0est = indirect0_haj,
                      ie0analytical = indirect0_analyticalvar_haj,
                      ie0boot = indirect0_bootvar_haj,
                      type = 'Hajek') %>% pivot_longer(ie0analytical:ie0boot)
ievartab2 = data.frame(gamma = gamma_list[2,],
                      #ie0est = indirect0_haj,
                      ie0analytical = indirect0_analyticalvar_ht,
                      ie0boot = indirect0_bootvar_ht,
                      type = 'HT') %>% pivot_longer(ie0analytical:ie0boot)
ievartab3 = data.frame(gamma = gamma_list[2,],
                       est = indirect0_haj,
                       ana = indirect0_analyticalvar_haj,
                       boot = indirect0_bootvar_haj,
                       type = 'HorvitzThompson') %>% 
  mutate(lb = est - 1.96*sqrt(boot),
         ub = est + 1.96*sqrt(boot)) 

ggplot(data = ievartab3, aes(x = gamma, y = est)) + 
  geom_point() + 
  geom_errorbar(aes(x = gamma, ymin = lb, ymax = ub))


ievartab = rbind(ievartab, ievartab2)

ggplot(ievartab,
       aes(x = gamma, y = value, group = name, color = name)) + 
  geom_line() + 
  facet_wrap(~type + name)

ievartab1 = data.frame(gamma = gamma_list[2,],
                       est = indirect1_haj,
                       boot = indirect1_bootvar_haj,
                       analytical = indirect1_analyticalvar_haj,
                       type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance')
ievartab2 = data.frame(gamma = gamma_list[2,],
                       est = indirect1_ht,
                       boot = indirect1_bootvar_ht,
                       analytical = indirect1_analyticalvar_ht,
                       type = 'Horvitz Thompson') %>% pivot_longer(boot:analytical, 'variance')
ievartab = rbind(ievartab1, ievartab2) %>%
  mutate(lb = est - 1.96*sqrt(value),
         ub = est + 1.96*sqrt(value),
         variance = ifelse(variance == 'boot', 'Bootstrapped Variance', 'Analytical Variance'))

ggplot(ievartab, aes(gamma, est)) + 
  geom_line() + facet_wrap(~type+variance) + 
  geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), fill="blue", alpha=0.25) + 
  labs(y = 'indirect effect',
       title = 'Indirect Effect Estimators and 95% Confidence Intervals',
       subtitle = '(treatment = 1)') + 
  theme(text = element_text(size = 14))



#OE VARIANCE COMPARISON
oetab = data.frame(gamma = gamma_list[2,],
                   oeanalytical = oe_analyticalvar_haj,
                   oeboot = oe_bootvar_haj, type = 'Hajek') %>% pivot_longer(oeanalytical:oeboot)
oetab2 = data.frame(gamma = gamma_list[2,],
                   oeanalytical = oe_analyticalvar_ht,
                   oeboot = oe_bootvar_ht, type = 'HorvitzThompson') %>% pivot_longer(oeanalytical:oeboot)
oetab = rbind(oetab, oetab2)
ggplot(oetab, aes(x = gamma, y = value, group = name, color = name)) + 
  geom_line() + 
  facet_wrap(~type + name)

oetab3 = data.frame(gamma = gamma_list[2,],
                    est = oe_haj,
                    lb = oe_haj - 1.96*sqrt(oe_bootvar_haj),
                    ub = oe_haj + 1.96*sqrt(oe_bootvar_haj))

ggplot(data = oetab3, aes(x = gamma, y = est)) + 
  geom_point() + 
  geom_errorbar(aes(x = gamma, ymin = lb, ymax = ub))

oevartab1 = data.frame(gamma = gamma_list[2,],
                       est = oe_haj,
                       boot = oe_bootvar_haj,
                       analytical = oe_analyticalvar_haj,
                       type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance')
oevartab2 = data.frame(gamma = gamma_list[2,],
                       est = oe_ht,
                       boot = oe_bootvar_ht,
                       analytical = oe_analyticalvar_ht,
                       type = 'Horvitz Thompson') %>% pivot_longer(boot:analytical, 'variance')
oevartab = rbind(oevartab1, oevartab2) %>%
  mutate(lb = est - 1.96*sqrt(value),
         ub = est + 1.96*sqrt(value),
         variance = ifelse(variance == 'boot', 'Bootstrapped Variance', 'Analytical Variance'))

ggplot(oevartab, aes(gamma, est)) + 
  geom_line() + facet_wrap(~type+variance) + 
  geom_errorbar(aes(x = gamma, ymin = lb, ymax = ub), fill="blue", alpha=0.25) + 
  labs(y = 'overall effect',
       title = 'Overall Effect Estimators and 95% Confidence Intervals') + 
  theme(text = element_text(size = 14))
 

comptab = rbind(oevartab %>% mutate(label = 'Overall Effect'),
                ievartab %>% mutate(label = 'Indirect Effect'),
                devartab %>% mutate(label = 'Direct Effect')) %>%
  filter(variance == 'Bootstrapped Variance')

data_hline <- data.frame(label = c('Indirect Effect', 'Direct Effect', 'Overall Effect'),  # Create data for lines
                         hline = c(0,3,0))
    

ggplot(comptab, aes(x = gamma, y = est)) + 
  geom_line() + facet_wrap(~type + label, scales = 'free_y') + 
  geom_ribbon(aes(x = gamma, ymin = lb, ymax =ub), fill = 'blue', alpha = 0.25) + 
  geom_hline(data = data_hline, aes(yintercept = hline, colour = 'red'), linetype = 'dashed') + 
  labs(y = 'causal estimates', title = 'Causal Estimates and 95% Confidence Intervals') + 
  theme(legend.position = 'none',
        text = element_text(size = 12))

