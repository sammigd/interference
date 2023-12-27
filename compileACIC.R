#compile big sim results
library(tidyverse)
setwd("/gpfs/gibbs/project/forastiere/sgd37/cai")
source('interference/compile_helper_funcs.R')
fig_loc = "~/project/cai/figures/ms_figs_bivar_one"

library(latex2exp)


#load("/gpfs/ysm/project/forastiere/sgd37/cai/figures/oct03wkspc.Rsave")

#save.image("/gpfs/ysm/project/forastiere/sgd37/cai/figures/nov13wkspc.Rsave")


#want separate analysis for each set of parameters
ngam = 100
gl = seq(from = -.2, to = .2, length.out = 33)
ngl = rep(0, 33)
gamma_list = rbind(rep(0, 99),
                   c(gl, ngl, ngl),
                   c(ngl, gl, ngl),
                   c(ngl, ngl, gl))
gamma_list = cbind(gamma_list, c(0,0,0,0))

beta_0 = .1
beta_1 = 3
#read in all of the simulations
#setwd("~/project/cai/florence_tfix5")
#setwd("~/project/cai/florence_tcheck")

#setwd("~/project/cai/powercalcs")
#setwd("~/project/cai/florence_g3_3")
#powercalcs results folder is 200 clusters, 15 members per cluster, epsilon ~ N(0,1)

setwd("~/project/cai/ms_stat_test_de_july3_univariate")


alliters = list.files()

#load(alliters[1])


#bigsims_homogtest_bigoe - folder with 200 sims of b3 = 3 to see if its a magnitude issue
#bigsims_threecovar - results for ACIC

#save.image("~/project/cai/bigsims_threecovar_full.Rsave")
#load("~/project/cai/bigsims_threecovar_full.Rsave")
#load(alliters[1])

#parameter sets 
table(str_remove(alliters, word(alliters, sep = '_')))
parmlist = unique(str_remove(alliters, word(alliters, sep = '_')))#[-1]
#parmlist = parmlist[7]

#will want to change this to looping over each of the parmlist elements i think
#iters = alliters[str_detect(alliters, parmlist[9])]

e.names = c('direct_ht', 'direct_analyticalvar_ht', 'direct_bootvar_ht',
            'direct_haj', 'direct_analyticalvar_haj', 'direct_bootvar_haj',
            'indirect0_ht', 'indirect0_analyticalvar_ht', 'indirect0_bootvar_ht',
            'indirect0_haj', 'indirect0_analyticalvar_haj', 'indirect0_bootvar_haj',
            'indirect1_ht', 'indirect1_analyticalvar_ht', 'indirect1_bootvar_ht',
            'indirect1_haj', 'indirect1_analyticalvar_haj', 'indirect1_bootvar_haj',
            'oe_haj', 'oe_analyticalvar_haj', 'oe_bootvar_haj',
            'oe_ht', 'oe_analyticalvar_ht', 'oe_bootvar_ht',
            'true_ie0', 'true_ie1', 'true_oe', 'true_de', 'true_y0', 'true_y1',
            'y0_ht', 'y1_ht', 'y0_haj', 'y1_haj', 'y0_ht_var', 'y1_ht_var', 'y0_ht_bootvar', 'y1_ht_bootvar', 'y0_ht_mcvar', 'y1_ht_mcvar', 
            'oe_anacoverage_ht', 'oe_anacoverage_haj', 'oe_bootcoverage_ht', 'oe_bootcoverage_haj',
            'de_anacoverage_ht', 'de_anacoverage_haj', 'de_bootcoverage_ht', 'de_bootcoverage_haj',
            'ie0_anacoverage_ht', 'ie0_anacoverage_haj', 'ie0_bootcoverage_ht', 'ie0_bootcoverage_haj',
            'ie1_anacoverage_ht', 'ie1_anacoverage_haj', 'ie1_bootcoverage_ht', 'ie1_bootcoverage_haj',
            'y0_anacoverage_ht', 'y1_anacoverage_ht', 'y_ht_covar', 'y0_bootcoverage_ht', 'y1_bootcoverage_ht', 'sum_wts0', 'sum_wts1'
)


#######################################
# POWER TESTING FOR FLAT OE ###########
#######################################
source('~/project/cai/interference/oe_stattest.R')

#get counts
table(sapply((str_split(alliters, '_', n = 2)), tail , 1))
sigresults = c()
empty = matrix(nrow = length(alliters), ncol = 3)
ind = 1

#for each set of parameters
for (pset in 1:length(parmlist)){
  iters = alliters[str_detect(alliters, parmlist[pset])]
  
  #for each simulation
  for (i in 1:length(iters)){
    load(iters[i])
    if(is.na(test$oe[1,1])){next}
    
    #run sig test for that simulation
    testresult = oe_sigtest(test$oe, test$oe_cov)
    #testresult_ie1 = oe_sigtest(test$ie1, test$ie1_cov)
    testresult_de = oe_sigtest(test$direct, test$direct_cov)
    
    sigresults = append(sigresults, testresult$accept)
    #maxlist = append(maxlist, testresult$diff)
    empty[ind,1] = parmlist[pset]
    empty[ind,2] = testresult$accept #is sig?
    empty[ind,3] = testresult$diff #max difference
    ind = ind + 1
  }
}

empty = data.frame(empty)
names(empty) = c('parms', 'TF', 'teststat')

summary_tab = data.frame(table(empty$parms, empty$TF)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  mutate(total = `FALSE` + `TRUE`,
         acc_rate =  `TRUE` / total,
         rej_rate = `FALSE` / total)
save(summary_tab, file = paste0(fig_loc, '/univar_power.RData') )


#########


bigbiastab = array(NA, dim = c(0, 30))
#compare bias for all the different parameters
for (parms in 1:length(parmlist)){
  print(parms)
  #parms = 1
  #make parms pretty
  ugly_parm = str_split(parmlist[parms], '_')
  pretty_parm = paste0('beta3 = ', ugly_parm[[1]][2], ', \nbeta4 = ', ugly_parm[[1]][3], ', \nconcor(X1,X2) = ', str_sub(ugly_parm[[1]][4], start = 1, end =-7 ))
  concordance =   str_sub(ugly_parm[[1]][4], start = 1, end =-7 )
  b3 = ugly_parm[[1]][2]
  b4 = ugly_parm[[1]][3]
  #if(!(b4 %in% c(0,1))){next}
  #if(!(concordance %in% 0)){next}
  estimates = get_means(parms, e.names)
  for (i in 1:length(e.names)){
    #print('newparms')
    assign(e.names[i], estimates[[i]])
  }
  biastab = data.frame(#g = gamma_list[2,]) %>% #
      g = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
      g1 = gamma_list[2,],
      g2 = gamma_list[3,],
      gamma_ind = c(rep('X1', floor(ngam/3)), rep('X2', floor(ngam/3)), rep('X3', floor(ngam/3)), 'X1')) %>% 
    mutate(#True parameter values
           true_oe = true_oe, true_ie1 = true_ie1, true_ie0 = true_ie0, 
           true_y0 = true_y0, true_y1 = true_y1, true_de = (true_y1 - true_y0),
           #estimated parameter values
           de_ht = direct_ht, de_haj = direct_haj,
           ie1_ht = indirect1_ht, ie0_ht = indirect0_ht, ie1_haj = indirect1_haj, ie0_haj = indirect0_haj,
           oe_ht = oe_ht, oe_haj = oe_haj,
           y1_ht = y1_ht, y1_haj=y1_haj, y0_ht = y0_ht, y0_haj = y0_haj,
           y0_ht_var = y0_ht_var, y1_ht_var = y1_ht_var,y_ht_covar = y_ht_covar,
           y0_ht_bootvar = y0_ht_bootvar, y1_ht_bootvar = y1_ht_bootvar,
           y0_ht_mcvar = y0_ht_mcvar, y1_ht_mcvar = y1_ht_mcvar, 
           #bias
           bias_de_ht = (direct_ht - beta_1), bias_de_haj = (direct_haj - beta_1),
           bias_oe_ht = (oe_ht - true_oe), bias_oe_haj = (oe_haj - true_oe),
           bias_ie0_ht = (indirect0_ht - true_ie0),  bias_ie1_ht = (indirect1_ht - true_ie1),
           bias_ie0_haj = (indirect0_haj - true_ie0), bias_ie1_haj = (indirect1_haj - true_ie1),
           bias_y0_ht = y0_ht - true_y0,
           bias_y1_ht = y1_ht - true_y1,
           bias_y0_haj = y0_haj - true_y0,
           bias_y1_haj = y1_haj - true_y1,
           #labels
           label = pretty_parm, concordance = concordance, b3 = b3, b4 = b4,
           #coverage
           oe_anacoverage_ht = oe_anacoverage_ht, oe_anacoverage_haj = oe_anacoverage_haj,
           oe_bootcoverage_ht = oe_bootcoverage_ht, oe_bootcoverage_haj = oe_bootcoverage_haj,
           de_anacoverage_ht = de_anacoverage_ht, de_anacoverage_haj = de_anacoverage_haj,
           de_bootcoverage_ht = de_bootcoverage_ht, de_bootcoverage_haj = de_bootcoverage_haj,
           ie0_anacoverage_ht = ie0_anacoverage_ht, ie0_anacoverage_haj = ie0_anacoverage_haj,
           ie0_bootcoverage_ht = ie0_bootcoverage_ht, ie0_bootcoverage_haj = ie0_bootcoverage_haj,
           ie1_anacoverage_ht = ie1_anacoverage_ht, ie1_anacoverage_haj = ie1_anacoverage_haj,
           ie1_bootcoverage_ht = ie1_bootcoverage_ht, ie1_bootcoverage_haj = ie1_bootcoverage_haj,
           y0_anacoverage_ht = y0_anacoverage_ht, y1_anacoverage_ht = y1_anacoverage_ht,
           y0_bootcoverage_ht = y0_bootcoverage_ht, y1_bootcoverage_ht = y1_bootcoverage_ht,
           #variances
           #OE HT
           oe_bootvar_ht = oe_bootvar_ht, 
           oe_analyticalvar_ht = oe_analyticalvar_ht,
           #OE HAJ
           oe_bootvar_haj = oe_bootvar_haj, 
           oe_analyticalvar_haj = oe_analyticalvar_haj,
           #DE HT
           de_analyticalvar_ht = direct_analyticalvar_ht, 
           de_bootvar_ht = direct_bootvar_ht,
           #DE HAJ
           de_analyticalvar_haj = direct_analyticalvar_haj,
           de_bootvar_haj = direct_bootvar_haj,
           #IE1 HT
           ie1_analyticalvar_ht = indirect1_analyticalvar_ht, 
           ie1_bootvar_ht = indirect1_bootvar_ht, 
           #IE1 HAJ
           ie1_analyticalvar_haj = indirect1_analyticalvar_haj, 
           ie1_bootvar_haj = indirect1_bootvar_haj,
           #IE0 HT
           ie0_analyticalvar_ht = indirect0_analyticalvar_ht, 
           ie0_bootvar_ht = indirect0_bootvar_ht, 
           #IE0 HAJ
           ie0_analyticalvar_haj = indirect0_analyticalvar_haj, 
           ie0_bootvar_haj = indirect0_bootvar_haj,
           #sum wts
           sum_wts0 = sum_wts0, sum_wts1 = sum_wts1
           )
  bigbiastab = rbind(bigbiastab, biastab)
} 

og = bigbiastab
#bigbiastab = og

#pick a for analytical variance, b for bootstrapped variance
bigbiastab = pick_var('a', bigbiastab)
bigbiastab = bigbiastab %>% filter(#gamma_ind =='X1', #removes a few dups
                                   (true_oe == 0 | g !='0'))

#######################################################
# FOR BIVARIATE INTERVENTION
#######################################################
bigbiastab$g
table(bigbiastab$gamma_ind)

ggplot(data = bigbiastab, aes(x = g1, y = g2, fill = oe_anacoverage_haj)) + 
  geom_tile() + 
  facet_wrap(~label, nrow = 2) + ggtitle('oe coverage for bivariate intervention')

ggplot(data = bigbiastab, aes(x = g1, y = oe_anacoverage_haj, col = g2)) + 
  geom_point() + 
  facet_wrap(~label) + ggtitle('oe coverage for bivariate intervention')

ggplot(data = bigbiastab, aes(x = g1, y = oe_bootcoverage_haj, col = g2)) + 
  geom_point() + 
  facet_wrap(~label) + ggtitle('oe coverage for bivariate intervention')


ggplot(data = bigbiastab, aes(x = g1, y = g2, fill = sum_wts0)) + 
  geom_tile() + 
  facet_wrap(~label)

ggplot(data = bigbiastab, aes(x = g1, y = g2, fill = oe_haj)) + 
  geom_tile() + 
  facet_wrap(~label, nrow = 2)


#######################################################
# CHECKING VARIANCES
#######################################################
#OE
vartab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, oe_analyticalvar_ht, oe_bootvar_ht, oe_analyticalvar_haj, oe_bootvar_haj) %>% 
  pivot_longer(oe_analyticalvar_ht:oe_bootvar_haj) 
ggplot(vartab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of Variances for OE') #+ ylim(0, 1E-5) #for zoomed in version
ggsave(paste0(fig_loc, '/mod_oe_vars.png'), width = 9, height = 9)


#DE
vartab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, de_analyticalvar_ht, de_analyticalvar_haj, de_bootvar_ht, de_bootvar_haj) %>% 
  pivot_longer(de_analyticalvar_ht:de_bootvar_haj) 
ggplot(vartab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of Variances for DE') #+ ylim(0, 1E-4) #for zoomed in version
ggsave(paste0(fig_loc, '/mod_de_vars.png'), width = 9, height = 9)


#ie
vartab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, ie1_analyticalvar_ht, ie1_analyticalvar_haj, ie1_bootvar_ht, ie1_bootvar_haj) %>% 
  pivot_longer(ie1_analyticalvar_ht:ie1_bootvar_haj) 
ggplot(vartab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of Variances for IE')
ggsave(paste0(fig_loc, '/mod_ie1_vars.png'), width = 9, height = 9)

vartab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, ie0_analyticalvar_ht, ie0_analyticalvar_haj, ie0_bootvar_ht, ie0_bootvar_haj) %>% 
  pivot_longer(ie0_analyticalvar_ht:ie0_bootvar_haj) 
ggplot(vartab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of Variances for IE')
ggsave(paste0(fig_loc, '/mod_ie0_vars.png'), width = 9, height = 9)




#yhats
vartab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, y0_ht_var, y1_ht_var, y0_ht_bootvar, y1_ht_bootvar, y0_ht_mcvar, y1_ht_mcvar) %>% 
  pivot_longer(y0_ht_var:y1_ht_mcvar) 
ggplot(vartab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of Variances for yhats')
ggsave(paste0(fig_loc, '/mod_yhat_ht_vars.png'), width = 9, height = 9)


#######################################################
# CHECKING COVERAGE
#######################################################
covtab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, oe_anacoverage_ht, oe_anacoverage_haj, 
                oe_bootcoverage_ht, oe_bootcoverage_haj) %>% 
  pivot_longer(oe_anacoverage_ht:oe_bootcoverage_haj) 
ggplot(covtab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Variance Coverage OE') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)
ggsave(paste0(fig_loc, '/mod_oe_cov1.png'), width = 9, height = 9)


ttt = bigbiastab %>% select(b3, b4, concordance, gamma_ind, g, true_ie1, ie1_ht, ie1_bootvar_ht)

ht_cov_tab = bigbiastab %>% 
  #filter(gamma_ind == 'X1') %>%
  dplyr::select(b3, b4, concordance, gamma_ind, g, 
                de_anacoverage_ht, ie0_anacoverage_ht, ie1_anacoverage_ht, oe_anacoverage_ht,
                de_anacoverage_haj, ie0_anacoverage_haj, ie1_anacoverage_haj, oe_anacoverage_haj) %>% 
  group_by(b3, b4, concordance, gamma_ind) %>% 
  summarise(DE_HT = round(mean(de_anacoverage_ht),2),
            IE0_HT = round(mean(ie0_anacoverage_ht),2),
            IE1_HT = round(mean(ie1_anacoverage_ht),2),
            OE_HT = round(mean(oe_anacoverage_ht),2),
            DE_HAJ = round(mean(de_anacoverage_haj),2),
            IE0_HAJ = round(mean(ie0_anacoverage_haj),2),
            IE1_HAJ = round(mean(ie1_anacoverage_haj),2),
            OE_HAJ = round(mean(oe_anacoverage_haj),2)) %>%
  mutate(DE = paste0(DE_HT, '  (', DE_HAJ, ')'),
         `IE(0)` =  paste0(IE0_HT, '  (', IE0_HAJ, ')'),
         `IE(1)` =  paste0(IE1_HT, '  (', IE1_HAJ,')'),
         OE = paste0(OE_HT, '  (', OE_HAJ, ')')) %>%
  dplyr::select(b3, b4, concordance, gamma_ind, DE, `IE(0)`, `IE(1)`, OE)
names(ht_cov_tab) =c('$beta3$', '$beta4$', 'Concordance(X1,X2)','Targeted Covariate', 'DE', 'IE(0)', 'IE(1)', 'OE')
#ht_cov_tab[,2:5] = round(ht_cov_tab[,2:5], 2)
#ht_cov_tab$Scenario = str_replace_all(ht_cov_tab$Scenario, 'beta3', '$beta_3$')
#ht_cov_tab$Scenario = str_replace_all(ht_cov_tab$Scenario, 'beta4', '$beta_4$')


library(xtable)
print(xtable(ht_cov_tab, type = "latex"), file = paste0(fig_loc, '/tabletest.tex'), sanitize.text.function = function(x){x}, include.rownames = F)



#IE0 coverage
covtab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, ie0_bootcoverage_ht, ie0_bootcoverage_haj, ie0_anacoverage_ht, ie0_anacoverage_haj) %>% 
  pivot_longer(ie0_bootcoverage_ht:ie0_anacoverage_haj) 
ggplot(covtab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Bootstrap Variance Coverage IE0') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)
ggsave(paste0(fig_loc, '/mod_ie0_cov1.png'), width = 9, height = 9)

#ie1
covtab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, ie1_anacoverage_ht, ie1_anacoverage_haj, ie1_bootcoverage_ht, ie1_bootcoverage_haj) %>% 
  pivot_longer(ie1_anacoverage_ht:ie1_bootcoverage_haj) 
ggplot(covtab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5)+ 
  facet_wrap(~label) + 
  ggtitle('Bootstrap Variance Coverage IE1') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)
ggsave(paste0(fig_loc, '/mod_ie1_cov1.png'), width = 9, height = 9)


#DE COV
covtab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj) %>% 
  pivot_longer(de_anacoverage_ht:de_bootcoverage_haj) 
ggplot(covtab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Variance Coverage DE') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)
ggsave(paste0(fig_loc, '/mod_de_cov2.png'), width = 9, height = 9)

#yhat coverage
covtab = bigbiastab %>% 
  dplyr::select(label, g, gamma_ind, y0_anacoverage_ht, y1_anacoverage_ht) %>% 
  pivot_longer(y0_anacoverage_ht:y1_anacoverage_ht) 
ggplot(covtab %>% filter(gamma_ind == 'X1'), aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Variance Coverage Yhat') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)
ggsave(paste0(fig_loc, '/mod_yhat_cov.png'), width = 9, height = 9)


#############################################
# POTENTIAL OUTCOMES
#############################################
y_bias_tab = bigbiastab %>%
  pivot_longer(cols = y1_ht:y0_haj, names_to = 'which_po', values_to = 'yhat')
ggplot(y_bias_tab %>% filter(gamma_ind == 'X1', concordance==0), aes(x =g, y = yhat, colour = which_po, shape = which_po)) + 
  geom_point() + 
  geom_line(aes(y = true_y0), colour = 'black') +
  geom_line(aes(y = true_y1), colour = 'black') + 
  facet_wrap(~label) #+ ylim(c(3,3.2))
ggsave(paste0(fig_loc, '/mod_yhats_est.png'), width = 7, height = 6)


y_bias_tab = bigbiastab %>%
  pivot_longer(cols = bias_y0_ht:bias_y1_haj, names_to = 'which_est', values_to = "bias") %>%
  filter(gamma_ind == 'X1', concordance==0)
ggplot(y_bias_tab,aes(x =g, y = bias, colour = which_est, shape = which_est)) + 
  geom_point() + 
  facet_wrap(~label)
ggsave(paste0(fig_loc, '/mod_yhats_bias.png'), width = 7, height = 6)




de_tab = bigbiastab %>%
  pivot_longer(cols = true_de:de_haj, names_to = 'which_de', values_to = 'de') 
  
ggplot(de_tab %>% filter(gamma_ind == 'X1', concordance==0), 
       aes(x =g, y = de, fill = which_de, colour = which_de)) +   
  geom_line() +
  #geom_line(aes(y = true_y0), colour = 'black') +
  #geom_line(aes(y = true_y1), colour = 'black') + 
  facet_wrap(~label + gamma_ind) + 
  #geom_ribbon(aes(x = g, ymin = de_ht_ana_lb, ymax = de_ht_ana_ub), alpha=0.2, fill = '#7CAE00', colour = NA) + 
  geom_ribbon(aes(x = g, ymin = de_haj_ana_lb, ymax = de_haj_ana_ub), alpha=0.2, fill = '#F8766D', colour = NA)
ggsave(paste0(fig_loc, '/mod_de.png'), width = 7, height = 6)


########################################################################
# FOR POWER CALCS
########################################################################
oe_bias_tab = bigbiastab %>%
  pivot_longer(cols = oe_ht:oe_haj, names_to = 'which_oe_est', values_to = 'oe_est') %>%
  pivot_longer(cols = c(oe_ht_ana_lb, oe_haj_ana_lb), names_to = "which_lb", values_to = 'lb') %>%
  pivot_longer(cols = c(oe_ht_ana_ub, oe_haj_ana_ub), names_to = "which_ub", values_to = 'ub') %>%
  filter(which_oe_est == 'oe_haj' & which_lb == 'oe_haj_ana_lb' |
           which_oe_est == 'oe_ht' & which_lb == 'oe_ht_ana_lb',
         which_oe_est == 'oe_haj' & which_ub == 'oe_haj_ana_ub' |
           which_oe_est == 'oe_ht' & which_ub == 'oe_ht_ana_ub')

#get oe true for power calcs
testing = oe_bias_tab %>% group_by(label, g) %>% summarise(true_oe = mean(true_oe))
testing = testing %>% group_by(label) %>% summarise(maxdiff = max(dist(true_oe)))


##########################################################
# EFFECT PLOTS FOR MS 
##########################################################
theme_set(theme_minimal()+
            theme(strip.text = element_text(size = 12),
                  text = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  legend.title = element_blank(),
                  legend.position = 'bottom',
                  panel.spacing = unit(1.4, "lines")))
for(i in c('DE', 'IE0', 'IE1', 'OE')){
  #i = 'OE'
  plot_df = make_ms_figtab(bigbiastab, effect = i)
  fig2df = plot_df$fig2df
  #ylabel = plot_df$ylab
  effect = plot_df$effect
  
  if(i == 'DE'){fig2df$true = fig2df$true_de; fig2df$est = fig2df$de_haj; ylab = 'Direct Effect'}
  if(i == 'IE0'){fig2df$true = fig2df$true_ie0; fig2df$est = fig2df$ie0_haj; ylab = 'Indirect Effect (0)'}
  if(i == 'IE1'){fig2df$true = fig2df$true_ie1; fig2df$est = fig2df$ie1_haj; ylab = 'Indirect Effect (1)'}
  if(i == 'OE'){fig2df$true = fig2df$true_oe; fig2df$est = fig2df$oe_haj; ylab = 'Overall Effect'}
  
  fig2df = fig2df %>% 
    filter(Conc == 'X1,X2 Uncorr' & gamma_ind == 'Target X1' | 
           Conc == 'X1,X2 Corr' & gamma_ind == 'Target X2' |
           Conc == 'X1,X2 Uncorr' & gamma_ind == 'Target X2') 
  
  #EFFECT PLOT
  ggplot(fig2df, aes(x = g, y = est))+#, #color = name, shape = name)) + 
    geom_point(colour = '#00BFC4') +
    geom_line(aes(y = true)) + 
    geom_ribbon(aes(x = g, ymin = lb, ymax = ub), 
                alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
    facet_grid(rows = vars(gamma_ind, Conc), cols = vars(Interference)) + 
    xlab(TeX(r'(\gamma)')) + 
    ylab(effect) #+ ylim(-.042,.042)
  ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_results.png'), width = 9, height = 6)
  
  #COVERAGE PLOT
  colname_boot = paste0(str_to_lower(i), '_bootcoverage_haj')
  colname_ana =paste0(str_to_lower(i), '_anacoverage_haj')
  
  covdf = fig2df %>% 
    pivot_longer(c(colname_boot,colname_ana)) %>% 
    mutate(name = case_when(name == colname_boot ~ 'bootstrap',
                            name == colname_ana ~ 'analytical'))
  
  ggplot(covdf, aes(x = g, y = value, colour = name)) + 
    geom_point(alpha = 0.5) + 
    facet_grid(rows = vars(gamma_ind, Conc), cols = vars(Interference)) + 
    #ggtitle('Coverage of 95% Confidence Intervals') +
    ylab('coverage') + 
    xlab(TeX(r'(\gamma)')) + 
    geom_hline(yintercept = 0.95)
  ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_coverage.png'), width = 9, height = 6)
  
  #BIAS PLOT
  colname_bias = paste0('bias_' ,str_to_lower(i), '_haj')
  fig2df$bias = fig2df[,colname_bias]
  ggplot(fig2df, aes(x = g, y = bias)) + 
    geom_point(alpha = 0.5) + 
    facet_grid(rows = vars(gamma_ind, Conc), cols = vars(Interference)) + 
    #ggtitle('Coverage of 95% Confidence Intervals') +
    ylab('bias') + 
    xlab(TeX(r'(\gamma)')) + 
    geom_hline(yintercept = 0)
  ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_bias.png'), width = 9, height = 6)
  
}


####################################################
# ALL SCENARIO PLOTS
####################################################
#OE
ggplot(bigbiastab %>% pivot_longer(c(true_oe, oe_ht)) ,aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 4, margin=margin(c(0,0,0,0))))+
  geom_ribbon(aes(x = g, ymin = oe_ht_ana_lb, ymax = oe_ht_ana_ub), alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
  facet_wrap(~gamma_ind+label, ncol = 6)
ggsave(paste0(fig_loc, '/mod_totalOE.png'), width = 9, height = 9)

#de 
ggplot(bigbiastab %>% pivot_longer(c(true_de, de_ht)) ,
       aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~gamma_ind + label, ncol = 4)
ggsave(paste0(fig_loc, '/mod_totalDE.png'), width = 9, height = 9)


#ie comp
ggplot(bigbiastab %>% pivot_longer(c(ie0_haj, true_ie0)) ,aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~gamma_ind + label, ncol = 4)
ggsave(paste0(fig_loc, '/mod_totalIE0.png'), width = 9, height = 9)

ggplot(bigbiastab %>% pivot_longer(c(ie1_haj, true_ie1)) %>%
         filter(gamma_ind == 'X1'),
       aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  geom_ribbon(aes(x = g, ymin = ie1_haj_ana_lb, ymax = ie1_haj_ana_ub), alpha=0.1) + 
  facet_wrap(~gamma_ind + label, ncol = 4)
ggsave(paste0(fig_loc, '/mod_totalIE1.png'), width = 9, height = 9)


#vary gammas for x1, x2, x3 faceted
devartab1 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = direct_haj,
                       boot = direct_bootvar_haj,
                       analytical = direct_analyticalvar_haj,
                       type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance') 
devartab2 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = direct_ht,
                       boot = direct_bootvar_ht,
                       analytical = direct_analyticalvar_ht,
                       type = 'Horvitz Thompson') %>% pivot_longer(boot:analytical, 'variance')
devartab = rbind(devartab1, devartab2) %>%
  mutate(lb = est - 1.96*sqrt(value),
         ub = est + 1.96*sqrt(value),
         variance = ifelse(variance == 'boot', 'Bootstrapped Variance', 'Analytical Variance'))

ggplot(devartab, aes(gamma, est)) + 
  geom_line() + facet_wrap(~type+variance+gamma_ind) + 
  geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), fill="blue", alpha=0.25) + 
  labs(y = 'direct effect',
       title = 'Direct Effect Estimators and 95% Confidence Intervals') + 
  theme(text = element_text(size = 14))

#indirect effect
ievartab1 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = indirect1_haj,
                       trueest = true_ie1,
                       boot = indirect1_bootvar_haj,
                       analytical = indirect1_analyticalvar_haj,
                       type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance')
ievartab2 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = indirect1_ht,
                       trueest = true_ie1,
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
  theme(text = element_text(size = 14)) + geom_line(aes(y = trueest), colour = 'red')




#overall effect
#OE VARIANCE COMPARISON
oevartab1 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = oe_haj,
                       trueest = true_oe,
                       boot = oe_bootvar_haj,
                       analytical = oe_analyticalvar_haj,
                       type = 'Hajek') %>% pivot_longer(boot:analytical, 'variance')
oevartab2 = data.frame(gamma = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
                       gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), NA),
                       est = oe_ht,
                       trueest = true_oe,
                       boot = oe_bootvar_ht,
                       analytical = oe_analyticalvar_ht,
                       type = 'Horvitz Thompson') %>% pivot_longer(boot:analytical, 'variance')
oevartab = rbind(oevartab1, oevartab2) %>%
  mutate(lb = est - 1.96*sqrt(value),
         ub = est + 1.96*sqrt(value),
         variance = ifelse(variance == 'boot', 'Bootstrapped Variance', 'Analytical Variance'))

ggplot(oevartab, aes(gamma, est)) + 
  geom_line() + facet_wrap(~type+variance+gamma_ind) + 
  geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), fill="blue", alpha=0.25) + 
  labs(y = 'overall effect',
       title = 'Overall Effect Estimators and 95% Confidence Intervals') + 
  theme(text = element_text(size = 14)) + 
  geom_line(aes(y=trueest), colour = 'red')




save.image(paste0(fig_loc, '/mod_wkspc.RSave'))

