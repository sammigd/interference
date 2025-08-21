#compile big sim results
#install.packages('tidyverse')
library(tidyverse)
library(viridis)
library(ggh4x)
library(reporter)
library(xtable)
library(here)
library(ggplot2)
library(latex2exp)


bivar = F
diffusion = F

setwd("/gpfs/gibbs/project/forastiere/sgd37/cai")
source('interference/compile_helper_funcs.R')

if(bivar){fig_loc = "~/project/cai/figures/ms_figs_bivar_two"}
#if(!bivar){fig_loc = "~/project/cai/figures/ms_figs_univar"}
if(!bivar){fig_loc = "~/project/cai/figures/ms_figs_univar_bigalpha"}
#if(!bivar){fig_loc = "~/project/cai/figures/ms_figs_univar_smallalpha"}



library(latex2exp)

theme_set(theme_minimal()+
            theme(strip.text = element_text(size = 10),
                  text = element_text(size = 11),
                  axis.text = element_text(size = 11),
                  legend.title = element_blank(),
                  legend.position = 'bottom',
                  panel.spacing = unit(.5, "lines"),
                  strip.background = element_rect(color = rgb(0,0,0,0))))


#load("/gpfs/ysm/project/forastiere/sgd37/cai/figures/oct03wkspc.Rsave")

#save.image("/gpfs/ysm/project/forastiere/sgd37/cai/figures/nov13wkspc.Rsave")


#univar heterogeneity
if(!bivar){
  ngam = 100
  gl = seq(from = -1.3, to = 1.3, length.out = 33)
  ngl = rep(0, 33)
  gamma_list = rbind(rep(0, 99),
                     c(gl, ngl, ngl),
                     c(ngl, gl, ngl),
                     c(ngl, ngl, gl))
  gamma_list = cbind(gamma_list, c(0,0,0,0))
  
  beta_0 = .1
  beta_1 = 3
}



#bivar heterogeneity
if(bivar){
  ngam = 64
  gl = c(seq(from = -1.3, to = 1.3, length.out = 7), 0)
  ggrid = expand.grid(gl, gl)
  gamma_list = t(cbind(0, ggrid, 0))
  
}

#read in all of the simulations
#setwd("~/project/cai/florence_tfix5")
#setwd("~/project/cai/florence_tcheck")

#setwd("~/project/cai/powercalcs")
#setwd("~/project/cai/florence_g3_3")
#powercalcs results folder is 200 clusters, 15 members per cluster, epsilon ~ N(0,1)

#ms stuff!!!
if(!bivar){setwd("~/project/cai/ms_univar_2025jul9_largealpha")}
#if(!bivar){setwd("~/project/cai/ms_univar_2025jul9_smallalpha")}

#if(!bivar){setwd("~/project/cai/ms_univar_2024may15")}

if(bivar){setwd("~/project/cai/ms_bivar_2024jun01")}



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
            'true_ie0', 'true_ie1', 'true_oe', 'true_de', 'true_y0', 'true_y1', 'true_y',
            'y0_ht', 'y1_ht', 'y0_haj', 'y1_haj', 'y_haj','y0_ht_var', 'y1_ht_var', 'y0_ht_bootvar', 'y1_ht_bootvar', 'y0_ht_mcvar', 'y1_ht_mcvar', 
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
sigresults_ie1 = c()
sigresults_ie0 = c()
sigresults_de = c()
empty = matrix(nrow = length(alliters), ncol = 3)
empty_other = matrix(nrow = length(alliters), ncol = 4)

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
    testresult_ie1 = oe_sigtest(test$indirect1, test$indirect_cov[,,2])
    testresult_ie0 = oe_sigtest(test$indirect0, test$indirect_cov[,,1])
    testresult_de = oe_sigtest(test$direct, test$direct_cov)
    
    sigresults = append(sigresults, testresult$accept)
    sigresults_ie1 = append(sigresults_ie1, testresult_ie1$accept)
    sigresults_ie0 = append(sigresults_ie0, testresult_ie0$accept)
    sigresults_de = append(sigresults_de, testresult_de$accept)
    
    
    #maxlist = append(maxlist, testresult$diff)
    empty[ind,1] = parmlist[pset]
    empty[ind,2] = testresult$accept #is sig?
    empty[ind,3] = testresult$diff #max difference
    
    empty_other[ind,1] = parmlist[pset]
    empty_other[ind,2] = testresult_ie0$accept
    empty_other[ind,3] = testresult_ie1$accept
    empty_other[ind,4] = testresult_de$accept
    
    ind = ind + 1
  }
  print(pset)
  print('is done')
  save(empty_other, file = paste0(fig_loc, '/bivar_power_ie_de.RData'))
  save(empty, file = paste0(fig_loc, '/bivar_power.RData'))

}

load(file = paste0(fig_loc, '/bivar_power.RData') )
load(file = paste0(fig_loc, '/bivar_power_ie_de.RData') )

empty_other = data.frame(empty_other)
names(empty_other) = c('parms', 'ie0_accept', 'ie1_accept', 'de_accept') 

de_power = empty_other %>% dplyr::select(parms, de_accept) %>%
  group_by(parms, de_accept) %>% summarise(Freq = n()) %>%
  pivot_wider(names_from = de_accept, values_from = Freq, values_fill = 0) %>%
  mutate(total = `FALSE` + `TRUE`,
         acc_rate =  `TRUE` / total,
         de_rej_rate = `FALSE` / total) %>% dplyr::select(parms, de_rej_rate)

ie0_power = empty_other %>% dplyr::select(parms, ie0_accept) %>%
  group_by(parms, ie0_accept) %>% summarise(Freq = n()) %>%
  pivot_wider(names_from = ie0_accept, values_from = Freq, values_fill = 0) %>%
  mutate(total = `FALSE` + `TRUE`,
         acc_rate =  `TRUE` / total,
         ie0_rej_rate = `FALSE` / total) %>% dplyr::select(parms, ie0_rej_rate)

ie1_power = empty_other %>% dplyr::select(parms, ie1_accept) %>%
  group_by(parms, ie1_accept) %>% summarise(Freq = n()) %>%
  pivot_wider(names_from = ie1_accept, values_from = Freq, values_fill = 0) %>%
  mutate(total = `FALSE` + `TRUE`,
         acc_rate =  `TRUE` / total,
         ie1_rej_rate = `FALSE` / total) %>% dplyr::select(parms, ie1_rej_rate)

other_power = merge(de_power, ie0_power, by = 'parms') %>% 
  merge(ie1_power, by = 'parms') %>%
  separate(parms, into = c('x', 'b3', 'b4', 'concordance', 'b5'), sep = '_') %>%
  mutate(b5 = substr(b5, 1,1)) %>%
  arrange(b4, b5, b3, concordance) %>%
  rename(de_level = de_rej_rate) %>%
  mutate(Ie0_Level = ifelse(b4==0 & b5==0, round(ie0_rej_rate,2), ' '),
         Ie0_Power = ifelse(b4!=0 | b5 !=0, round(ie0_rej_rate,2), ' '),
         Ie1_Level = ifelse(b4==0 & b5==0, round(ie1_rej_rate,2), ' '),
         Ie1_Power = ifelse(b4!=0 | b5 !=0, round(ie1_rej_rate,2), ' ')) %>%
  dplyr::select(b3, b4, b5, concordance, de_level, Ie0_Level, Ie0_Power, Ie1_Level, Ie1_Power)

print(xtable(other_power, type = "latex"), file = "../power_table_ie_de.tex", include.rownames = F)


empty = data.frame(empty)
names(empty) = c('parms', 'TF', 'teststat')
str(empty)
teststattable = empty %>% 
  mutate(teststat = as.numeric(teststat)) %>%
  group_by(parms) %>% summarise(teststat = mean(teststat, na.rm = T))

summary_tab = data.frame(table(empty$parms, empty$TF)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  mutate(total = `FALSE` + `TRUE`,
         acc_rate =  `TRUE` / total,
         rej_rate = `FALSE` / total)# %>%
  #separate(Var1, into = c('x', 'b3', 'b4', 'concordance', 'b5'), sep = '_') %>%
  #mutate(b5 = substr(b5, 1,1))
save(summary_tab, file = paste0(fig_loc, '/univar_power.RData') )

power_tab = teststattable %>% 
  merge(summary_tab, by.x = 'parms', by.y = 'Var1')  %>%
  separate(parms, into = c('x', 'b3', 'b4', 'concordance', 'b5'), sep = '_') %>%
  mutate(b5 = substr(b5, 1,1)) %>% dplyr::select(-x, -`TRUE`, -`FALSE`, -total) %>%
  arrange(b4, b5, b3, concordance) %>%
  mutate(Level = ifelse(b4==0 & b5==0, round(rej_rate,2), ' '),
         Power = ifelse(b4!=0 | b5 !=0, round(rej_rate,2), ' ')) %>%
  dplyr::select(b3, b4, b5, concordance, teststat, Level, Power)
  
print(xtable(power_tab, type = "latex"), file = "power_table.tex", include.rownames = F)



#########
#if(!bivar){beta5 = 0; b5 = 0}

bigbiastab = array(NA, dim = c(0, 30))
#compare bias for all the different parameters
for (parms in 1:length(parmlist)){
  print(parms)
  #make parms pretty
  ugly_parm = str_split(parmlist[parms], '_')
  pretty_parm = paste0('beta3 = ', ugly_parm[[1]][2], ', \nbeta4 = ', ugly_parm[[1]][3], ', \nconcor(X1,X2) = ', str_sub(ugly_parm[[1]][4], start = 1, end =-7 ))
  concordance = ugly_parm[[1]][4]
  b3 = ugly_parm[[1]][2]
  b4 = ugly_parm[[1]][3]
  b5 = str_sub(ugly_parm[[1]][5], 1, 1)
  #if(!bivar){b5 = 0}
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
           label = pretty_parm, concordance = concordance, b3 = b3, b4 = b4, b5 = b5,
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

hist(bigbiastab$true_de)

og = bigbiastab 
#bigbiastab = og

#pick a for analytical variance, b for bootstrapped variance
bigbiastab = pick_var('a', bigbiastab) 

if(!bivar){
  #bigbiastab = bigbiastab %>% 
  #  filter(gamma_ind == 'X1' & (true_oe == 0 | g !='0')) %>%
  #  filter(gamma_ind != 'X3')
  bigbiastab = bigbiastab %>% filter(#gamma_ind =='X1', #removes a few dups
    (true_oe == 0 | g !='0')) %>%
    filter(b5 == 0)
}
#
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
# EFFECT / COVERAGE / BIAS PLOTS FOR MS 
##########################################################

for(i in c('DE', 'IE0', 'IE1', 'OE')){
  #i = 'OE'
  plot_df = make_ms_figtab(bigbiastab, effect = i, diffusion = F)
  fig2df = plot_df$fig2df
  
  #ylabel = plot_df$ylab
  effect = plot_df$effect
  
  if(i == 'DE'){fig2df$true = fig2df$true_de; fig2df$est = fig2df$de_haj; ylab = 'Direct Effect'}
  if(i == 'IE0'){fig2df$true = fig2df$true_ie0; fig2df$est = fig2df$ie0_haj; ylab = 'Indirect Effect (0)'}
  if(i == 'IE1'){fig2df$true = fig2df$true_ie1; fig2df$est = fig2df$ie1_haj; ylab = 'Indirect Effect (1)'}
  if(i == 'OE'){fig2df$true = fig2df$true_oe; fig2df$est = fig2df$oe_haj; ylab = 'Overall Effect'}
 
  a = 'Cor(X' %p% supsc('(1)') %p%',X' %p% supsc('(2)') %p% ')=0'
  b = 'Cor(X' %p% supsc('(1)') %p%',X' %p% supsc('(2)') %p% ')>0'
  
  if(!bivar){
    fig2df = fig2df %>% 
      filter(Conc == a & gamma_ind == 'Target X' %p% supsc('(1)') | 
             Conc == b & gamma_ind == 'Target X' %p% supsc('(2)') |
             Conc == a & gamma_ind == 'Target X' %p% supsc('(2)')) 
  }
  
  fig2df = fig2df %>%
    mutate(IntIn1 = 'Interference Variable: X' %p% supsc('(1)'),#'Interference Variable: X1', #expression(paste(X^(1))) is right expression doesnt work 
           IntIn2 = 'Interference Variable: X' %p% supsc('(2)'))
  
  
  #EFFECT PLOT
  if(bivar){
    dd = fig2df %>%
      filter(!(Interference == 'No Interference\n (β₃=0; β₄=0; β₅=0)' & 
               InterferenceX2 =='Heterogeneous Interference\n (β₅=1)'),
             !(Interference == 'No Interference\n (β₃=0; β₄=0; β₅=0)' & 
                 Conc =='Cor(X⁽¹⁾,X⁽²⁾)>0')) %>%
      mutate(InterferenceX2 = case_when(InterferenceX2 == 'No Interference\n (β₅=0)' ~ 'Homogeneous Interference\n (β₅=0)',
                                        InterferenceX2 == 'Heterogeneous Interference\n (β₅=1)' ~ 'Heterogeneous Interference\n (β₅=1)'),
             InterferenceX2 = factor(InterferenceX2, levels = c('Homogeneous Interference\n (β₅=0)','Heterogeneous Interference\n (β₅=1)' )),
             IntIn1 = ifelse(Interference == 'No Interference\n (β₃=0; β₄=0; β₅=0)', 'No Interference', IntIn1),
             IntIn1 = fct_rev(factor(IntIn1)))#,
             #IntIn1 = ifelse(Interference == 'No Interference\n (β₃=0; β₄=0)', '', IntIn1))
    ggplot(data = dd, aes(x = g1, y = g2, fill = est)) + 
      geom_raster() + 
      ggh4x::facet_nested(rows = vars(IntIn2, InterferenceX2, Conc), 
                          cols = vars(IntIn1, Interference), 
                          nest_line = element_line(), 
                          solo_line = T,
                          render_empty = F) + 
      xlab(expression(paste(gamma[1]))) +
      ylab(expression(paste(gamma[2]))) +
      theme(legend.text = element_text(angle=-90, hjust = -.05)) + 
      scale_fill_viridis() + 
      theme(legend.title = element_text(size = 12)) + 
      labs(fill = effect)
    #ggsave(paste0(fig_loc, '/', 'new_OE_labels.png'), width = 9, height = 9)
    
    ggsave(paste0(fig_loc, '/', effect, '_bivar_modelsim_results.png'), width = 9, height = 9)
  }else{ #univar
    ggplot(fig2df, aes(x = g, y = est))+#, #color = name, shape = name)) + 
      geom_point(colour = '#00BFC4') +
      geom_line(aes(y = true)) + 
      geom_ribbon(aes(x = g, ymin = lb, ymax = ub), 
                  alpha = .2, fill = '#00BFC4', colour = rgb(0,0,0,0)) +
      facet_nested(rows = vars(gamma_ind, Conc), 
                   cols = vars(IntIn1, Interference), 
                   nest_line = element_line(), 
                   solo_line = T,
                   render_empty = F) + 
      xlab(TeX(r'(\gamma)')) + 
      ylab(ylab) #+ ylim(-.042,.042)
    ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_results.png'), width = 9, height = 6)
  }
  
  #COVERAGE PLOT
  colname_boot = paste0(str_to_lower(i), '_bootcoverage_haj')
  colname_ana =paste0(str_to_lower(i), '_anacoverage_haj')
  
  if(bivar){
    covdf = dd %>% 
      pivot_longer(c(colname_boot,colname_ana)) %>% 
      mutate(name = case_when(name == colname_boot ~ 'bootstrap',
                              name == colname_ana ~ 'analytical'))
    
    ggplot(data = covdf %>% filter(name == 'analytical'), aes(x = g1, y = g2, fill = value)) + 
      geom_raster() + 
      facet_nested(rows = vars(IntIn2, InterferenceX2, Conc), 
                   cols = vars(IntIn1, Interference), 
                   nest_line = element_line(), 
                   solo_line = T,
                   render_empty = F) + 
      xlab(expression(paste(gamma[1]))) +
      ylab(expression(paste(gamma[2]))) +
      scale_fill_viridis() + 
      theme(legend.title = element_text(size = 12),
            legend.text = element_text(angle=-90, hjust = -.05)) +
      labs(fill = '95% CI coverage')
    
    ggsave(paste0(fig_loc, '/', effect, '_bivar_modelsim_coverage.png'), width = 9, height = 9)
  }else{
    covdf = fig2df %>% 
      pivot_longer(c(colname_boot,colname_ana)) %>% 
      mutate(name = case_when(name == colname_boot ~ 'bootstrap',
                              name == colname_ana ~ 'analytical'))
    
    ggplot(covdf, aes(x = g, y = value, colour = name)) + 
      geom_point(alpha = 0.5) + 
      facet_nested(rows = vars(gamma_ind, Conc), cols = vars(IntIn1, Interference), nest_line = element_line(), solo_line = T) + 
      ylab('95% CI Coverage') + 
      xlab(TeX(r'(\gamma)')) + 
      geom_hline(yintercept = 0.95)
    ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_coverage.png'), width = 9, height = 6)
  }
  
  #BIAS PLOT

  
  if(bivar){
    colname_bias = paste0('bias_' ,str_to_lower(i), '_haj')
    dd$bias = dd[,colname_bias]
    ggplot(data = dd, aes(x = g1, y = g2, fill = bias)) + 
      geom_raster() + 
      facet_nested(rows = vars(IntIn2, InterferenceX2, Conc), 
                   cols = vars(IntIn1, Interference), nest_line = element_line(), solo_line = T, render_empty = F) + 
      xlab(expression(paste(gamma[1]))) +
      ylab(expression(paste(gamma[2]))) +
      scale_fill_viridis() + 
      theme(legend.title = element_text(size = 12),
            legend.text = element_text(angle=-90, hjust = -.05)) + 
      labs(fill = paste0(effect, ' Bias'))
    ggsave(paste0(fig_loc, '/', effect, '_bivar_modelsim_bias.png'), width = 9, height = 9)
    
  }else{
    colname_bias = paste0('bias_' ,str_to_lower(i), '_haj')
    fig2df$bias = fig2df[,colname_bias]
    ggplot(fig2df, aes(x = g, y = bias)) + 
      geom_point(alpha = 0.5) + 
      facet_nested(rows = vars(gamma_ind, Conc), 
                   cols = vars(IntIn1, Interference), 
                   nest_line = element_line(), 
                   solo_line = T, 
                   render_empty = F) + 
      ylab('Bias') + 
      xlab(TeX(r'(\gamma)')) + 
      geom_hline(yintercept = 0)
    ggsave(paste0(fig_loc, '/', effect, '_univar_modelsim_bias.png'), width = 9, height = 6)
  }

}

#bias tables for appendix
##bivariate
bt = bigbiastab %>%
  filter(gamma_ind == 'X1') %>%
  group_by(b3,b4,b5, concordance) %>%
  summarise(DE = round(mean(bias_de_haj),6),
            IE0 = round(mean(bias_ie0_haj),6),
            IE1 = round(mean(bias_ie1_haj),6),
            OE = round(mean(bias_oe_haj),6)) %>%
  #mutate(gamma_ind = ifelse(gamma_ind == 'X1', 'Target X' %p% supsc('(1)'), 'Target X' %p% supsc('(2)'))) %>%
  arrange(b3, b4, b5, concordance)

print(xtable(bt, type = "latex", digits = 6, align = "|c|c|c|c|c|c|c|c|c|"), file = here("bias_table.tex"), include.rownames = F)



##univariate
uni_bt = bigbiastab %>%
  filter(gamma_ind != 'X3') %>%
  mutate(concordance = ifelse(concordance == '0.Rsave', "0", "0.65"),
         b5 = '0') %>%
  group_by(b3,b4,b5,concordance, gamma_ind) %>%
  summarise(`DE bias` = round(mean(bias_de_haj),6),
            `IE(0) bias` = round(mean(bias_ie0_haj),6),
            `IE(1) bias` = round(mean(bias_ie1_haj),6),
            `OE bias` = round(mean(bias_oe_haj),6)) %>%
  #mutate(gamma_ind = ifelse(gamma_ind == 'X1', 'Target X' %p% supsc('(1)'), 'Target X' %p% supsc('(2)'))) %>%
  arrange(b3, b4, b5, concordance)

print(xtable(uni_bt, type = "latex", digits = 6, align = "|c|c|c|c|c|c|c|c|c|c|"), file = here("uni_bias_table.tex"), include.rownames = F)


save.image(paste0(fig_loc, '/mod_wkspc.RSave'))

