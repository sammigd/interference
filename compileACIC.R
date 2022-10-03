#compile big sim results
library(tidyverse)
setwd("/gpfs/ysm/project/forastiere/sgd37/cai")
source('scripts/compile_helper_funcs.R')
#library(geepack)

#load("/gpfs/ysm/project/forastiere/sgd37/cai/figures/sep29wkspc.Rsave")

#save.image("/gpfs/ysm/project/forastiere/sgd37/cai/figures/oct03wkspc.Rsave")


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
setwd("~/project/cai/florence_n5_v2")
#setwd("~/project/cai/florence")
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
#parmlist = parmlist[1:18]

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
            'true_ie0', 'true_ie1', 'true_oe', 'y0_ht', 'y1_ht', 'y0_haj', 'y1_haj',
            'oe_anacoverage_ht', 'oe_anacoverage_haj', 'oe_bootcoverage_ht', 'oe_bootcoverage_haj',
            'de_anacoverage_ht', 'de_anacoverage_haj', 'de_bootcoverage_ht', 'de_bootcoverage_haj',
            'oe_truecov_ht', 'oe_truevar_ht', 'oe_mcvar_ht', 'oe_mccoverage_ht') 

#######################################
# POWER TESTING FOR FLAT OE ###########
#######################################
source('~/project/cai/scripts/oe_stattest.R')

#get counts
table(sapply((str_split(alliters, '_', n = 2)), tail , 1))

empty = matrix(nrow = length(alliters), ncol = 3)
ind = 1
for (pset in 1:length(parmlist)){
  iters = alliters[str_detect(alliters, parmlist[pset])]
  for (i in 1:length(iters)){
    load(iters[i])
    if(is.na(test$oe[1,1])){next}
    testresult = oe_sigtest(test$oe, test$oe_cov)
    sigresults = append(sigresults, testresult$accept)
    maxlist = append(maxlist, testresult$diff)
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
  estimates = get_means(parms, e.names)
  for (i in 1:length(e.names)){
    assign(e.names[i], estimates[[i]])
  }
  biastab = data.frame(#g = gamma_list[2,]) %>% #
      g = c(rep(gamma_list[2,1:floor(ngam/3)], 3), 0),
      gamma_ind = c(rep('X1', 33), rep('X2', 33), rep('X3', 33), 'X1')) %>% 
    mutate(true_oe = true_oe, true_ie1 = true_ie1, true_ie0 = true_ie0, true_de = beta_1,
           de_ht = direct_ht, de_haj = direct_haj,
           ie1_ht = indirect1_ht, ie0_ht = indirect0_ht, ie1_haj = indirect1_haj, ie0_haj = indirect0_haj,
           oe_ht = oe_ht, oe_haj = oe_haj,
           bias_de_ht = (direct_ht - beta_1), bias_de_haj = (direct_haj - beta_1),
           bias_oe_ht = (oe_ht - true_oe), bias_oe_haj = (oe_haj - true_oe),
           bias_ie0_ht = (indirect0_ht - true_ie0),  bias_ie1_ht = (indirect1_ht - true_ie1),
           bias_ie0_haj = (indirect0_haj - true_ie0), bias_ie1_haj = (indirect1_haj - true_ie1),
           label = pretty_parm, concordance = concordance,
           y1_ht = y1_ht, y1_haj=y1_haj, y0_ht = y0_ht, y0_haj = y0_haj,
           oe_anacoverage_ht = oe_anacoverage_ht, oe_anacoverage_haj = oe_anacoverage_haj,
           oe_bootcoverage_ht = oe_bootcoverage_ht, oe_bootcoverage_haj = oe_bootcoverage_haj,
           de_anacoverage_ht = de_anacoverage_ht, de_anacoverage_haj = de_anacoverage_haj,
           de_bootcoverage_ht = de_bootcoverage_ht, de_bootcoverage_haj = de_bootcoverage_haj,
           oe_truecov_ht = oe_truecov_ht, oe_truevar_ht = oe_truevar_ht, oe_bootvar_ht = oe_bootvar_ht, oe_analyticalvar_ht = oe_analyticalvar_ht,
           oe_mcvar_ht = oe_mcvar_ht, oe_mccoverage_ht = oe_mccoverage_ht
           )
  bigbiastab = rbind(bigbiastab, biastab)
} 
bigbiastab2 = bigbiastab %>% pivot_longer(cols = bias_de_ht:bias_ie0_haj, names_to = 'which_estimator', values_to = 'bias')

#bias investigation
bias_tab = bigbiastab %>% filter(gamma_ind == "X1") 
ggplot(bias_tab, aes(x = g, y = bias_oe_ht)) + 
  geom_point(alpha = .4) + 
  facet_wrap(~label) + 
  geom_hline(yintercept = 0) + 
  xlab('gamma') + 
  ylab('oe bias (haj)') + 
  title('Bias In Hajek Overall Effect')


boxplot(bias_tab$bias_oe_ht ~ bias_tab$g)
abline(h=0, col = 'red')
boxplot(bigbiastab$bias_oe_haj ~ bigbiastab$g)
abline(h=0, col = 'red')


#bias for each set of params
ggplot(bigbiastab2, aes(x = g, y = bias)) + 
  geom_point() + 
  facet_wrap(~which_estimator + label, ncol = length(parmlist))

og = bigbiastab #flo 400
og2 = bigbiastab #flo 200
bigbiastab = og2

#pick a for analytical variance, b for bootstrapped variance
bigbiastab = pick_var('b', bigbiastab)
bigbiastab = bigbiastab %>% filter(gamma_ind =='X1', (true_oe == 0 | g !='0'), concordance ==0)
test = bigbiastab %>% filter(g == '0', gamma_ind =='X1')

table(bigbiastab$g)

#coverage'
cov_tab = bigbiastab %>%
  select(label, g, gamma_ind, oe_truecov_ht) %>% distinct()
table(bigbiastab$gamma_ind)

pdf("/gpfs/ysm/project/forastiere/sgd37/cai/figures/covtab_400clusters.png")  
boxplot(cov_tab$oe_truecov_ht ~ cov_tab$g, xlab = 'gamma', ylab = 'oe coverage', main = 'oe ht true variance coverage with 200 clusters')
abline(h = 0.95 , col = 'red')
abline(h = 0.5 , col = 'pink')
dev.off()


#comparing variances
vartab = bigbiastab %>% 
  dplyr::select(label, g, oe_analyticalvar_ht, oe_bootvar_ht, oe_mcvar_ht) %>% 
  pivot_longer(oe_analyticalvar_ht:oe_mcvar_ht) 
ggplot(vartab, aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of True, Bootstrapped, and Analytical Variance for OE, HT
          \n beta4=0 -> trueoe=0 and n_i = 5 and no rnorm') #+ ylim(0, 1E-4) #for zoomed in version

#comparing coverages
covtab = bigbiastab %>% 
  dplyr::select(label, g, oe_anacoverage_ht, oe_bootcoverage_ht, oe_mccoverage_ht) %>% 
  pivot_longer(oe_anacoverage_ht:oe_mccoverage_ht) 
ggplot(covtab, aes(x = g, y = value, colour = name)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~label) + 
  ggtitle('Comparison of True, Bootstrapped, and MC Coverage for OE, HT
          \n beta4=0 -> trueoe=0 and n_i = 5 and no rnorm') +
  ylab('coverage') + 
  geom_hline(yintercept = 0.95)



#potential outcomes
y_bias_tab = bigbiastab %>%
  pivot_longer(cols = y1_ht:y0_haj, names_to = 'which_po', values_to = 'yhat')
ggplot(y_bias_tab %>% filter(gamma_ind == 'X1', concordance==0), aes(x =g, y = yhat, colour = which_po)) + 
  geom_point() + 
  facet_wrap(~label)

#direct effect bias
de_bias_tab = bigbiastab %>%
  pivot_longer(cols = de_ht:de_haj, names_to = 'which_de_est', values_to = 'de_est') %>%
  pivot_longer(cols = c(de_ht_ana_lb, de_haj_ana_lb), names_to = "which_lb", values_to = 'lb') %>%
  pivot_longer(cols = c(de_ht_ana_ub, de_haj_ana_ub), names_to = "which_ub", values_to = 'ub') %>%
  filter(which_de_est == 'de_haj' & which_lb == 'de_haj_ana_lb' |
           which_de_est == 'de_ht' & which_lb == 'de_ht_ana_lb',
         which_de_est == 'de_haj' & which_ub == 'de_haj_ana_ub' |
           which_de_est == 'de_ht' & which_ub == 'de_ht_ana_ub')

ggplot(de_bias_tab %>% filter(gamma_ind == 'X1'), aes(x = g, y = bias_de_ht, fill = which_de_est)) + 
  geom_line(aes(colour = which_de_est)) + 
  geom_ribbon(aes(x = g, ymin = lb, ymax = ub), alpha=0.25) + 
  facet_wrap(~label, ncol = length(parmlist)) + 
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') + 
  theme(strip.text = element_text(size = 4)) +
  labs(title = 'Direct Effect Bias - HT and Haj')

#indirect effect bias
ie0_bias_tab = bigbiastab %>%
  pivot_longer(cols = ie0_ht:ie0_haj, names_to = 'which_ie0_est', values_to = 'ie0_est') %>%
  pivot_longer(cols = c(ie0_ht_ana_lb, ie0_haj_ana_lb), names_to = "which_lb", values_to = 'lb') %>%
  pivot_longer(cols = c(ie0_ht_ana_ub, ie0_haj_ana_ub), names_to = "which_ub", values_to = 'ub') %>%
  filter(which_ie0_est == 'ie0_haj' & which_lb == 'ie0_haj_ana_lb' |
           which_ie0_est == 'ie0_ht' & which_lb == 'ie0_ht_ana_lb',
         which_ie0_est == 'ie0_haj' & which_ub == 'ie0_haj_ana_ub' |
           which_ie0_est == 'ie0_ht' & which_ub == 'ie0_ht_ana_ub')

ggplot(ie0_bias_tab, aes(x = g, y = bias_ie0_ht, fill = which_ie0_est)) + 
  geom_point(aes(colour = which_ie0_est)) + 
  geom_ribbon(aes(x = g, ymin = lb, ymax = ub), alpha=0.25) + 
  facet_wrap(~label, ncol = length(parmlist)/2) + 
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') + 
  theme(strip.text = element_text(size = 7)) + 
  labs(title = 'Indirect Effect (Trt = 0) Bias - HT and Haj')

ie1_bias_tab = bigbiastab %>%
  pivot_longer(cols = ie1_ht:ie1_haj, names_to = 'which_ie1_est', values_to = 'ie1_est') %>%
  pivot_longer(cols = c(ie1_ht_ana_lb, ie1_haj_ana_lb), names_to = "which_lb", values_to = 'lb') %>%
  pivot_longer(cols = c(ie1_ht_ana_ub, ie1_haj_ana_ub), names_to = "which_ub", values_to = 'ub') %>%
  filter(which_ie1_est == 'ie1_haj' & which_lb == 'ie1_haj_ana_lb' |
           which_ie1_est == 'ie1_ht' & which_lb == 'ie1_ht_ana_lb',
         which_ie1_est == 'ie1_haj' & which_ub == 'ie1_haj_ana_ub' |
           which_ie1_est == 'ie1_ht' & which_ub == 'ie1_ht_ana_ub')

ggplot(ie1_bias_tab, aes(x = g, y = bias_ie1_ht, fill = which_ie1_est)) + 
  geom_line(aes(colour = which_ie1_est)) + 
  geom_ribbon(aes(x = g, ymin = lb, ymax = ub), alpha=0.25) + 
  facet_wrap(~label, ncol = length(parmlist)) + 
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') + 
  theme(strip.text = element_text(size = 7)) + 
  labs(title = 'Indirect Effect (Trt = 1) Bias - HT and Haj')


#OVERALL EFFECT
oe_bias_tab = bigbiastab %>%
  pivot_longer(cols = oe_ht:oe_haj, names_to = 'which_oe_est', values_to = 'oe_est') %>%
  pivot_longer(cols = c(oe_ht_ana_lb, oe_haj_ana_lb), names_to = "which_lb", values_to = 'lb') %>%
  pivot_longer(cols = c(oe_ht_ana_ub, oe_haj_ana_ub), names_to = "which_ub", values_to = 'ub') %>%
  filter(which_oe_est == 'oe_haj' & which_lb == 'oe_haj_ana_lb' |
           which_oe_est == 'oe_ht' & which_lb == 'oe_ht_ana_lb',
         which_oe_est == 'oe_haj' & which_ub == 'oe_haj_ana_ub' |
           which_oe_est == 'oe_ht' & which_ub == 'oe_ht_ana_ub')

ggplot(oe_bias_tab %>% filter(gamma_ind == 'X1'), aes(x = g, y = bias_oe_ht, fill = which_oe_est)) + 
  geom_line(aes(colour = which_oe_est)) + 
  geom_ribbon(aes(x = g, ymin = lb, ymax = ub), alpha=0.25) + 
  facet_wrap(~label+gamma_ind, ncol = length(parmlist)/2) + 
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') + 
  theme(strip.text = element_text(size = 5))+
  labs(title = 'Overall Effect Bias - HT and Haj')

#get oe true for power calcs
testing = oe_bias_tab %>% group_by(label, g) %>% summarise(true_oe = mean(true_oe))
testing = testing %>% group_by(label) %>% summarise(maxdiff = max(dist(true_oe)))

#Does the "truth" make sense
ggplot(bigbiastab %>% pivot_longer(true_oe:true_de),aes(x = g, y = value)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~name+label, ncol = 18)

#Do the estimates make sense
ggplot(bigbiastab %>% pivot_longer(c(ie0_haj, ie1_haj)) %>%filter(gamma_ind == 'X1') ,aes(x = g, y = value)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~name+label, ncol = 9)+
  ylim(-.1, .1)

library(ggtext)

#wht does confidence interval look the same for every scenario?



#oe comp
fig2df = bigbiastab %>% pivot_longer(c(true_oe, oe_haj)) %>%
  #filter(gamma_ind == 'X1') %>%
  filter(str_detect(label, 'beta3 = 0'), concordance != "0.65") %>%
  mutate(Interference = case_when(str_detect(label, 'beta4 = 0,') ~ 'No Interference',
                                  str_detect(label, 'beta4 = 0.5') ~ 'Moderate Interference',
                                  str_detect(label, 'beta4 = 1') ~ 'Strong Interference')) %>%
  filter(!is.na(Interference)) %>%
  mutate(Conc = case_when(concordance == '0' ~ 'X1,X2 Uncorrelated',
                          concordance == '0.8' ~ 'X1,X2 Highly Correlated')) %>%
  mutate(Interference = factor(Interference, levels = c('No Interference',
                                                        'Moderate Interference',
                                                        'Strong Interference'))) %>%
  mutate(Conc = factor(Conc, levels = c('X1,X2 Uncorrelated', 'X1,X2 Highly Correlated')),
         name = ifelse(name == 'oe_haj', 'Hajek OE', 'True OE'),
         name = factor(name, levels = c( 'True OE','Hajek OE')))

#install.packages('latex2exp')
library(latex2exp)

ggplot(fig2df %>% filter(gamma_ind == 'X1'),aes(x = g, y = value, colour = name)) + 
  theme_minimal() +
  geom_line() +   
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = 'bottom')+
  geom_ribbon(aes(x = g, ymin = oe_haj_ana_lb, ymax = oe_haj_ana_ub), 
              alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
  facet_grid(rows = vars(Conc), cols = vars(Interference)) + 
  xlab(TeX(r'(\gamma)'))+ 
  ylab('Overall Effect') + ylim(-.042,.042)
#ggsave('x1intervention.png', width = 7, height = 5)
#library(Cairo)
#cairo_pdf(file = "~/project/cai/ggplot-greek.pdf", width = 8, height = 5)
## ggplot object created here
#dev.off()

ggplot(fig2df %>% filter(gamma_ind == 'X2'),aes(x = g, y = value, colour = name)) + 
  theme_minimal() +
  geom_line() +   
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = 'bottom')+
  geom_ribbon(aes(x = g, ymin = oe_haj_ana_lb, ymax = oe_haj_ana_ub), 
              alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
  facet_grid(rows = vars(Conc), cols = vars(Interference)) + 
  xlab('gamma') + ylab('Overall Effect') + ylim(-.042,.042)

#1100 x 800

#testing ci
ggplot(bigbiastab %>% filter(gamma_ind == 'X1'), aes(x = g, y = oe_haj)) + 
  geom_line() + 
  geom_errorbar(aes(x = g, ymin = oe_haj_ana_lb, ymax = oe_haj_ana_lb), colour = 'red') + 
  facet_wrap(~label)


ggplot(bigbiastab %>% pivot_longer(c(true_oe, oe_ht)) ,aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 4, margin=margin(c(0,0,0,0))))+
  geom_ribbon(aes(x = g, ymin = oe_ht_ana_lb, ymax = oe_ht_ana_ub), alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
  facet_wrap(~gamma_ind+label, ncol = 18)

#de comp
ggplot(bigbiastab %>% pivot_longer(c(true_de, de_ht)) ,
       aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~gamma_ind + label, ncol = 18)

#ie comp
ggplot(bigbiastab %>% pivot_longer(c(ie0_haj, true_ie0)) ,aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~gamma_ind + label, ncol = 9)
ggplot(bigbiastab %>% pivot_longer(c(ie1_haj, true_ie1)) ,aes(x = g, y = value, colour = name)) + 
  geom_line() +   
  theme(strip.text = element_text(size = 7, margin=margin(c(0,0,0,0))))+
  facet_wrap(~gamma_ind + label, ncol = 9)

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



