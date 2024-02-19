library(stringr)
library(latex2exp)
library(tidyverse)
library(reporter)

#load the helper functions
source('/home/sgd37/project/cai/interference/compile_helper_funcs.R')

#load the truth
load('simtruthwkspc.Rsave')

#set the output folder
fig_loc = "~/project/cai/figures/nov16_laura"
figloc2 = 'nov16_laura'

#set the folder where the simulation results are
setwd("~/project/cai/ms_diffusion_2023nov17_wtruth")

#load in the diffusion truth simulation results
#load('~/project/cai/Y0_truth_sim_df.Rsave')


#get the parameters
alliters = list.files()
table(str_remove(alliters, word(alliters, sep = '_')))
parmlist = unique(str_remove(alliters, word(alliters, sep = '_')))#[-7]


gl = seq(from = -.2, to = .2, length.out = 33)
ngl = rep(0, 33)
gamma_list = rbind(rep(0, 33),
                   c(ngl, gl),
                   c(gl, ngl))
gamma_list = cbind(gamma_list, c(0,0, 0))
gamma_numer = gamma_list
ngam = ncol(gamma_numer)
diffusion = T

b5 = NULL
beta5 = NULL
beta4 = NULL

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

#load('/gpfs/gibbs/project/forastiere/sgd37/cai/simtruthwkspc.Rsave')

bigbiastab = array(NA, dim = c(0, 30))
#compare bias for all the different parameters
for (parms in 1:length(parmlist)){
  print(paste('parms=',parms))
  #parms = 1
  #make parms pretty
  ugly_parm = str_split(parmlist[parms], '_')
  pretty_parm = paste0('conc = ', ugly_parm[[1]][2], ', \np(diff) = ', str_sub(ugly_parm[[1]][3], start = 1, end =-7 ))
  #print(pretty_parm)
  concordance = ugly_parm[[1]][2]
  pdiff = str_sub(ugly_parm[[1]][3], start = 1, end =-7 )
  truth_use = Y0_truth_bar[Y0_truth_bar$concordance == concordance & Y0_truth_bar$pdiff == pdiff,] %>% 
    mutate(g = as.numeric(g)) %>%
    arrange(g) %>%
    arrange(desc(gamma_ind))
  estimates = get_means(parms, e.names, beta4 =Inf, diffusion = T, truth_use = truth_use)
  for (i in 1:length(e.names)){
    assign(e.names[i], estimates[[i]])
  }
  biastab = data.frame(#g = gamma_list[2,], #) %>% #
    g = c(rep(gamma_list[3,1:floor(ngam/2)], 2), 0),
    gamma_ind = c(rep('X2', 33), rep('X1', 33), 'X1')) %>% 
    mutate(#True parameter values
           true_oe = true_oe, 
           #true_ie1 = true_ie1, true_ie0 = true_ie0, 
           #true_y0 = true_y0, true_y1 = true_y1, true_de = 0,#(true_y1 - true_y0),
           #estimated parameter values
           de_ht = direct_ht, de_haj = direct_haj,
           ie1_ht = indirect1_ht, ie0_ht = indirect0_ht, ie1_haj = indirect1_haj, ie0_haj = indirect0_haj,
           oe_ht = oe_ht, oe_haj = oe_haj,
           y1_ht = y1_ht, y1_haj=y1_haj, y0_ht = y0_ht, y0_haj = y0_haj,
           #bias
           #bias_de_ht = (direct_ht - beta_1), bias_de_haj = (direct_haj - beta_1),
           #bias_oe_ht = (oe_ht - true_oe), bias_oe_haj = (oe_haj - true_oe),
           #bias_ie0_ht = (indirect0_ht - true_ie0),  bias_ie1_ht = (indirect1_ht - true_ie1),
           #bias_ie0_haj = (indirect0_haj - true_ie0), bias_ie1_haj = (indirect1_haj - true_ie1),
           #labels
           label = pretty_parm, concordance = concordance, pdiff = pdiff,
           #coverage
           oe_anacoverage_ht = oe_anacoverage_ht, oe_anacoverage_haj = oe_anacoverage_haj,
           oe_bootcoverage_ht = oe_bootcoverage_ht, oe_bootcoverage_haj = oe_bootcoverage_haj,
           de_anacoverage_ht = de_anacoverage_ht, de_anacoverage_haj = de_anacoverage_haj,
           de_bootcoverage_ht = de_bootcoverage_ht, de_bootcoverage_haj = de_bootcoverage_haj,
           ie0_anacoverage_ht = ie0_anacoverage_ht, ie0_anacoverage_haj = ie0_anacoverage_haj,
           ie0_bootcoverage_ht = ie0_bootcoverage_ht, ie0_bootcoverage_haj = ie0_bootcoverage_haj,
           ie1_anacoverage_ht = ie1_anacoverage_ht, ie1_anacoverage_haj = ie1_anacoverage_haj,
           ie1_bootcoverage_ht = ie1_bootcoverage_ht, ie1_bootcoverage_haj = ie1_bootcoverage_haj,
           #variances
           oe_bootvar_ht = oe_bootvar_ht, oe_analyticalvar_ht = oe_analyticalvar_ht,
           oe_bootvar_haj = oe_bootvar_haj, oe_analyticalvar_haj = oe_analyticalvar_haj,
           de_analyticalvar_ht = direct_analyticalvar_ht, de_analyticalvar_haj = direct_analyticalvar_haj,
           de_bootvar_ht = direct_bootvar_ht, de_bootvar_haj = direct_bootvar_haj,
           #sum of weights
           sum_wts0 = sum_wts0, sum_wts1 = sum_wts1
    )
  bigbiastab = rbind(bigbiastab, biastab)
} 

bigbiastab = pick_var('a', bigbiastab) #this makes CIs

#bigbiastab = merge(bigbiastab, Y0_truth_bar, by = c('concordance', 'pdiff', 'g', 'gamma_ind'), all.x = T)
#bigbiastab = merge(bigbiastab, true_oe_df, by = c('concordance', 'pdiff', 'g', 'gamma_ind'), all.x = T)

names(bigbiastab)

#plot y0
ggplot(bigbiastab %>% filter(concordance == 0),
       aes(x = g, y = y0_haj, colour = concordance, fill = concordance)) +
  geom_point() + 
  geom_line(aes(y = true_y0), colour = 'black') + 
  #geom_ribbon(aes(x = g, ymin = de_haj_ana_lb, ymax = de_haj_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(gamma_ind), rows = vars(pdiff))
ggsave(here('figures', 'figloc2', 'diffusiony0_truth.png'), width = 6, height = 6)

#plot y1
ggplot(bigbiastab %>% filter(concordance == 0),
       aes(x = g, y = y1_haj, colour = concordance, fill = concordance)) +
  geom_point() + 
  geom_line(aes(y = true_y1), colour = 'black') + 
  #geom_ribbon(aes(x = g, ymin = de_haj_ana_lb, ymax = de_haj_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(gamma_ind), rows = vars(pdiff))
ggsave(here('figures', 'figloc2', 'diffusiony1_truth.png'), width = 6, height = 6)

#plot direct effect
ggplot(bigbiastab,
       aes(x = g, y = de_ht, colour = concordance, fill = concordance)) +
  geom_point() + 
  geom_line(aes(y = true_de), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = de_ht_ana_lb, ymax = de_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(gamma_ind), rows = vars(pdiff))
ggsave(paste0(fig_loc, '/diffusionde.png'), width = 6, height = 6)


ggplot(bigbiastab %>% 
         select(g, gamma_ind, pdiff, concordance, de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj) %>%
         pivot_longer(de_anacoverage_ht:de_bootcoverage_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  facet_wrap(~concordance+gamma_ind+pdiff)
ggsave(paste0(fig_loc, '/decoverage.png'), width = 6, height = 6)

plot_df = make_ms_figtab(bigbiastab, effect = 'OE', diffusion = T)
fig2df = plot_df$fig2df
#plot overall effect
ggplot(fig2df,
       aes(x = g, y = oe_haj)) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.spacing = unit(1.4, "lines")) +
  geom_point(colour = '#00BFC4') +
  geom_line(aes(y = true_oe)) + 
  #geom_line(true_oe_df, aes(x = g, y = true_oe)) + 
  geom_ribbon(aes(x = g, ymin = lb, ymax = ub), 
              alpha = 0.2,fill = '#00BFC4', colour = rgb(0,0,0,0)) + 
  facet_grid(rows = vars(gamma_ind, Conc), cols = vars(pdiff_long)) + 
  xlab(TeX(r'(\gamma)')) + ylab('Overall Effect (Hajek Estimator)')

ggplot(true_oe_df, aes(x = g, y = true_oe)) + 
  facet_grid(rows = vars(gamma_ind, concordance), cols = vars(pdiff)) +
  geom_line() + 
  geom_point(data = fig2df, aes(x = g, y = oe_haj))




ggsave(paste0(fig_loc, '/diffusionoe.png'), width = 8, height = 6)

#oe coverage
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, oe_anacoverage_ht, oe_anacoverage_haj, oe_bootcoverage_ht, oe_bootcoverage_haj) %>%
         pivot_longer(oe_anacoverage_ht:oe_bootcoverage_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  geom_line(aes(y = 0.95),colour = 'black') + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusion_oe_coverage.png'), width = 6, height = 6)


#plot ie
ggplot(bigbiastab, 
       aes(x = g, y = ie0_haj, colour = concordance, fill = concordance)) +
  geom_point() + 
  geom_line(aes(y = true_ie0), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = ie0_ht_ana_lb, ymax = ie0_ht_ana_ub), alpha=0.25) + 
  facet_grid(rows = vars(gamma_ind), cols = vars(pdiff))
ggsave(paste0(fig_loc, '/diffusion_ie0.png'), width = 6, height = 6)


ggplot(bigbiastab,
       aes(x = g, y = ie1_haj, colour = concordance, fill = concordance)) +
  geom_point() + 
  geom_line(aes(y = true_ie1), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = ie1_ht_ana_lb, ymax = ie1_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(gamma_ind), rows = vars(pdiff))
ggsave(paste0(fig_loc, '/diffusion_ie1.png'), width = 6, height = 6)


#ie0 coverage
ggplot(bigbiastab %>% 
         select(g, gamma_ind, pdiff,concordance, ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj) %>%
         pivot_longer(ie0_anacoverage_ht:ie0_bootcoverage_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  facet_grid(cols = vars(concordance, pdiff), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusion_ie0_cov.png'), width = 6, height = 6)

save.image(paste0(fig_loc, 'diffusionwkspc.RSave'))
