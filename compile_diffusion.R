source('/gpfs/ysm/project/forastiere/sgd37/cai/scripts/compile_helper_funcs.R')

fig_loc = "~/project/cai/figures/nov23meeting"

setwd("~/project/cai/diffusion1")
alliters = list.files()
parmlist = unique(str_sub(alliters, -7, -7))
beta4 = Inf
ngam = ncol(gamma_numer)

e.names = c('direct_ht', 'direct_analyticalvar_ht', 'direct_bootvar_ht',
            'direct_haj', 'direct_analyticalvar_haj', 'direct_bootvar_haj',
            'indirect0_ht', 'indirect0_analyticalvar_ht', 'indirect0_bootvar_ht',
            'indirect0_haj', 'indirect0_analyticalvar_haj', 'indirect0_bootvar_haj',
            'indirect1_ht', 'indirect1_analyticalvar_ht', 'indirect1_bootvar_ht',
            'indirect1_haj', 'indirect1_analyticalvar_haj', 'indirect1_bootvar_haj',
            'oe_haj', 'oe_analyticalvar_haj', 'oe_bootvar_haj',
            'oe_ht', 'oe_analyticalvar_ht', 'oe_bootvar_ht',
            'true_ie0', 'true_ie1', 'true_oe', 'true_de', 'true_y0', 'true_y1',
            'y0_ht', 'y1_ht', 'y0_haj', 'y1_haj',
            'oe_anacoverage_ht', 'oe_anacoverage_haj', 'oe_bootcoverage_ht', 'oe_bootcoverage_haj',
            'de_anacoverage_ht', 'de_anacoverage_haj', 'de_bootcoverage_ht', 'de_bootcoverage_haj', 
            'ie0_anacoverage_ht', 'ie0_anacoverage_haj', 'ie0_bootcoverage_ht', 'ie0_bootcoverage_haj',
            'oe_mcvar_ht') 

bigbiastab = array(NA, dim = c(0, 30))
#compare bias for all the different parameters
for (parms in 1:length(parmlist)){
  print(parms)
  #parms = 1
  #make parms pretty
  ugly_parm = str_split(parmlist[parms], '_')
  pretty_parm = ugly_parm
  #pretty_parm = paste0('beta3 = ', ugly_parm[[1]][2], ', \nbeta4 = ', ugly_parm[[1]][3], ', \nconcor(X1,X2) = ', str_sub(ugly_parm[[1]][4], start = 1, end =-7 ))
  #concordance =   str_sub(ugly_parm[[1]][4], start = 1, end =-7 )
  concordance = parmlist[parms]
  estimates = get_means(parms, e.names, beta4 =Inf, diffusion = T)
  for (i in 1:length(e.names)){
    assign(e.names[i], estimates[[i]])
  }
  biastab = data.frame(#g = gamma_list[2,]) %>% #
    g = c(rep(gamma_list[3,1:floor(ngam/2)], 2), 0),
    gamma_ind = c(rep('X2', 33), rep('X1', 33), 'X1')) %>% 
    mutate(true_oe = true_oe, true_ie1 = true_ie1, true_ie0 = true_ie0, true_de = true_de,
           true_y0 = true_y0, true_y1 = true_y1,
           de_ht = direct_ht, de_haj = direct_haj,
           ie1_ht = indirect1_ht, ie0_ht = indirect0_ht, ie1_haj = indirect1_haj, ie0_haj = indirect0_haj,
           oe_ht = oe_ht, oe_haj = oe_haj,
           bias_de_ht = (direct_ht - true_de), bias_de_haj = (direct_haj - true_de),
           bias_oe_ht = (oe_ht - true_oe), bias_oe_haj = (oe_haj - true_oe),
           bias_ie0_ht = (indirect0_ht - true_ie0),  bias_ie1_ht = (indirect1_ht - true_ie1),
           bias_ie0_haj = (indirect0_haj - true_ie0), bias_ie1_haj = (indirect1_haj - true_ie1),
           label = pretty_parm, concordance = concordance,
           y1_ht = y1_ht, y1_haj=y1_haj, y0_ht = y0_ht, y0_haj = y0_haj,
           oe_anacoverage_ht = oe_anacoverage_ht, oe_anacoverage_haj = oe_anacoverage_haj,
           oe_bootcoverage_ht = oe_bootcoverage_ht, oe_bootcoverage_haj = oe_bootcoverage_haj,
           de_anacoverage_ht = de_anacoverage_ht, de_anacoverage_haj = de_anacoverage_haj,
           de_bootcoverage_ht = de_bootcoverage_ht, de_bootcoverage_haj = de_bootcoverage_haj,
           ie0_anacoverage_ht = ie0_anacoverage_ht, ie0_anacoverage_haj = ie0_anacoverage_haj, 
           ie0_bootcoverage_ht = ie0_bootcoverage_ht, ie0_bootcoverage_haj = ie0_bootcoverage_haj,

           oe_bootvar_ht = oe_bootvar_ht, oe_analyticalvar_ht = oe_analyticalvar_ht,
           oe_mcvar_ht = oe_mcvar_ht
    )
  bigbiastab = rbind(bigbiastab, biastab)
} 

bigbiastab = pick_var('b', bigbiastab)

names(bigbiastab)

bigbiastab = bigbiastab %>%
  mutate(concordance = case_when(concordance == 0 ~ 'concordance 0.5',
                                 concordance == 5 ~ 'concordance 0.65'))

#plot y0
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, y0_ht, y0_haj, true_y0) %>%
         pivot_longer(y0_ht:y0_haj), 
       aes(x = g, y = value, colour = name, fill = name)) +
  geom_point() + 
  geom_line(aes(y = true_y0), colour = 'black') + 
  #geom_ribbon(aes(x = g, ymin = de_ht_ana_lb, ymax = de_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusiony0.png'), width = 6, height = 6)

#plot y1
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, y1_ht, y1_haj, true_y1) %>%
         pivot_longer(y1_ht:y1_haj), 
       aes(x = g, y = value, colour = name, fill = name)) +
  geom_point() + 
  geom_line(aes(y = true_y1), colour = 'black') + 
  #geom_ribbon(aes(x = g, ymin = de_ht_ana_lb, ymax = de_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusiony1.png'), width = 6, height = 6)

#plot direct effect
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, de_ht, de_haj, true_de, de_ht_ana_ub, de_ht_ana_lb) %>%
         pivot_longer(de_ht:de_haj), 
       aes(x = g, y = value, colour = name, fill = name)) +
  geom_point() + 
  geom_line(aes(y = true_de), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = de_ht_ana_lb, ymax = de_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusionde.png'), width = 6, height = 6)


ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, de_anacoverage_ht, de_anacoverage_haj, de_bootcoverage_ht, de_bootcoverage_haj) %>%
         pivot_longer(de_anacoverage_ht:de_bootcoverage_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  facet_wrap(~concordance+gamma_ind)


#plot overall effect
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, oe_ht, oe_haj, true_oe, oe_ht_ana_lb, oe_ht_ana_ub) %>%
         pivot_longer(oe_ht:oe_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  geom_line(aes(y = true_oe),colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = oe_ht_ana_lb, ymax = oe_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusionoe.png'), width = 6, height = 6)

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
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, ie0_ht, ie0_haj, true_ie0, ie0_ht_ana_lb, ie0_ht_ana_ub) %>%
         pivot_longer(ie0_ht:ie0_haj), 
       aes(x = g, y = value, colour = name, fill = name)) +
  geom_point() + 
  geom_line(aes(y = true_ie0), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = ie0_ht_ana_lb, ymax = ie0_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusion_ie0.png'), width = 6, height = 6)


ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, ie1_ht, ie1_haj, true_ie1, ie1_ht_ana_lb, ie1_ht_ana_ub) %>%
         pivot_longer(ie1_ht:ie1_haj), 
       aes(x = g, y = value, colour = name, fill = name)) +
  geom_point() + 
  geom_line(aes(y = true_ie1), colour = 'black') + 
  geom_ribbon(aes(x = g, ymin = ie1_ht_ana_lb, ymax = ie1_ht_ana_ub), alpha=0.25) + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusion_ie1.png'), width = 6, height = 6)


#ie0 coverage
ggplot(bigbiastab %>% 
         select(g, gamma_ind, concordance, ie0_anacoverage_ht, ie0_anacoverage_haj, ie0_bootcoverage_ht, ie0_bootcoverage_haj) %>%
         pivot_longer(ie0_anacoverage_ht:ie0_bootcoverage_haj), 
       aes(x = g, y = value, fill = name, colour = name)) +
  geom_point() + 
  facet_grid(cols = vars(concordance), rows = vars(gamma_ind))
ggsave(paste0(fig_loc, '/diffusion_ie0_cov.png'), width = 6, height = 6)

save.image(paste0(fig_loc, 'diffusionwkspc.RSave'))
