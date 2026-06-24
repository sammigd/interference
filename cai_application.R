library(here)
here()
library(tidyverse)
library(abind)
library(latex2exp)
library(viridis)
library(rlist)

repo_loc = "~/project_pi_lf474/sgd37/"

source(paste0(repo_loc, 'interference/load_clean_cai.R')) 

library(Interference, lib.loc = paste0(repo_loc, '/pckgs'))


source(paste0(repo_loc, 'interference/', 'analysis_scripts/GroupIPW_sd v2.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/denominator_sd.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/helper_functs.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/parallel_bootvar_function test.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/ypop_sd.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/hajek_adj.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/de_sd.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/ie_sd.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/CalcHeterogenousTrueIE.R')) #flag loading issues 2026/4/4
source(paste0(repo_loc, 'interference/', 'analysis_scripts/oe_sd.R'))
source(paste0(repo_loc, 'interference/', 'analysis_scripts/oe_stattest.R'))

revisions = T
seed = 812

fig_loc = 'ms_figs_cai' #for main publication figs
fig_loc_revisions = 'revision_figs_application'
if(revisions == T){fig_loc = fig_loc_revisions}

#check correlation between degree and median rank of neighbors degree
df = df %>% group_by(cluster) %>% mutate(dgr_rank = rank(dgr))
ids = V(cai_graph)
df$mean_dgr_rank = NA
for (i in ids){
  #get vertex
  i_use = ids[i]
  
  #get neighbors of vertex
  my_neigh = neighbors(cai_graph, v = i_use)
  i_use = as.numeric(names(i_use))
  
  #get avg rank of neighbors of vertex
  df$mean_dgr_rank[df$id == i_use] <-  mean(df$dgr_rank[df$id %in% as.numeric(names(my_neigh))])
  
}

cor(df$dgr, df$mean_dgr_rank, use = 'complete.obs')


plot(df$dgr, df$mean_dgr_rank)

#check which clusters are big / small
clustersize = df %>% group_by(cluster) %>% summarise(size = n(), num_trt = sum(trt)) %>%
  mutate(alpha = num_trt / size)
plot(clustersize$size, clustersize$alpha)


df = merge(df, clustersize, by = 'cluster') %>% 
  group_by(cluster) %>% 
  mutate(avg_degree = mean(norm_degree))

#compare distributution of centrality inn big v small clusters
df$cluster_cat = ifelse(df$cluster %in% clustersize$cluster[clustersize$size <= 80], '<=80', '>80')

df = df %>% mutate(cluster_cat = case_when(size < 30 ~ '<30',
                                           size >=30 & size < 50 ~ '30-49',
                                           size >=50 & size < 80 ~ '50-79',
                                           size>=80 ~ '80+'))

ggplot(df, aes(x = norm_btwn, color = cluster_cat)) + geom_density() + 
  #ggtitle('Density of Norm Betweenness in Large and Small Clusters') + 
  theme_bw() +
  labs(x = 'Normalized Betweenness', y = 'Density', colour = 'Village Size')
#ggsave(here('figures', fig_loc, paste0('comparing_normbtwn_split.png')), height = 4, width = 6)

ggplot(df, aes(x = norm_degree, color = cluster_cat)) + geom_density() + 
  #ggtitle('Density of Norm Degree in Large and Small Clusters') + 
  theme_bw() +
  labs(x = 'Normalized Degree', y = 'Density', colour = 'Village Size')
#ggsave(here('figures', fig_loc, paste0('comparing_normdgr_split.png')), height = 4, width = 6)


ggplot(df, aes(x = dgr, color = cluster_cat)) + geom_density() + ggtitle('Density of Degree in Large and Small Clusters')
#ggsave(here('figures', fig_loc, paste0('comparing_dgr_split.png')), height = 5, width = 5)


ggplot(df, aes(x = n.peers, fill = cluster_cat)) + geom_bar() + ggtitle('Density of n.peers in Large and Small Clusters')
#ggsave(here('figures', fig_loc, paste0('comparing_npeers_split.png')), height = 5, width = 5)


#check if it works if I impute missing values from group mean or something like that
names(df)
#22 - general_trust
#26 - understanding
#17 - ricearea_2010

table(is.na(df$general_trust))
table(is.na(df$understanding))
table(is.na(df$ricearea_2010))
table(is.na(df$risk_averse))
table(is.na(df$disaster_prob))
table(is.na(df$disaster_loss))
table(is.na(df$educ_good))



df$understanding_imp = df$understanding
df$understanding_imp[is.na(df$understanding_imp)] <- mean(df$understanding, na.rm = T)

df$ricearea_2010_imp = df$ricearea_2010
df$ricearea_2010_imp[is.na(df$ricearea_2010_imp)] <- mean(df$ricearea_2010, na.rm = T)

df$educ_good_imp = df$educ_good
df$educ_good_imp[is.na(df$educ_good_imp)] <- mean(df$educ_good, na.rm = T)

df$age_imp = df$age
df$age_imp[is.na(df$age_imp)] <- mean(df$age, na.rm = T)

df$male_imp = df$male
df$male_imp[is.na(df$male_imp)] <- mean(df$male, na.rm = T)

df$lit_imp = df$literacy
df$lit_imp[is.na(df$literacy)] <- mean(df$literacy, na.rm = T)



names(df)
bivar = F

olddf = df
#df = olddf
trt_col_list = c(#55, #educ_good
  54, #rice area 2010
  21, #disaster prob
  #28, #risk averse
  34,#, #dgr
  35#, #btwn
  #53, #understanding
  #56 #age
)
#trt_col_list = list(c(34,35), c(34,54)) #bivariate
if(bivar){trt_col_list = list(c(34,54), c(34,35))}
#trt_col_list = list(c(34,54))
de_results = list()
ie0_results = list()
ie1_results = list()
oe_results = list()
univar_pvals = matrix(NA, nrow = 0, ncol = 6)
bivar_pvals = matrix(NA, nrow = 0, ncol = 8)
bivar_oe_results = list()
de_univar_pvals = matrix(NA, nrow = 0, ncol = 6)
ie0_univar_pvals = matrix(NA, nrow = 0, ncol = 6)
ie1_univar_pvals = matrix(NA, nrow = 0, ncol = 6)

univar_gamma_Rs = c(25, 55, 100, 199)
bivar_gamma_Rs = c(2, 4, 7, 7*2, 7*3, 7*4)

j_range = 1:3
#if(revisions == T){j_range = 3}

for(gamma_R in univar_gamma_Rs){
  counter = 0
  for(j in j_range){ #1:3
    #j = 1
    df = olddf
    if (j == 1){
      ww = '_over80'
      df = df[df$cluster %in% clustersize$cluster[clustersize$size > 80],]
    }
    
    if(j == 2){
      ww = '_le80'
      df = df[df$cluster %in% clustersize$cluster[clustersize$size <= 80],]
    }
    
    if(j==3){
      ww = '_all'
      df = olddf
    }
    
    relabel_key = data.frame(old_cl = unique(df$cluster), new_cl = 1:length(unique(df$cluster)))
    df = merge(df, relabel_key, by.x = 'cluster', by.y = 'old_cl')
    df$neigh = df$new_cl
    
    for (i in trt_col_list){
      counter = counter+1
      #i = 21
      #i = c(34, 35)
      estimand <- '1'
      cov_cols = c(i) #cov_cols = c(34,54)
      
      trt_col = 41
      out_col = 42
      #ps_with_re = FALSE
      #numerator_with_re = FALSE
      out_name = 'trt'
      if(!bivar){target_name = paste0(names(df)[i], ww)} #target_name = 'dgr_btwn'
      if(bivar){target_name = paste0('dgr_ricearea',ww)}
      
      clean_name = case_when(target_name %in% c('btwn_all', 'btwn_le80', 'btwn_over80') ~ 'Betweenness',
                             target_name %in% c('dgr_all' , 'dgr_le80' , 'dgr_over80') ~ 'Degree',
                             target_name %in% c('disaster_prob_all' , 'disaster_prob_le80' , 'disaster_prob_over80') ~ 'Future Disaster Probability',
                             target_name %in% c('ricearea_2010_imp_all' , 'ricearea_2010_imp_le80' , 'ricearea_2010_imp_over80') ~ 'Rice Production Area')
      
      phi_hat <- list(coefs = c(0, 0), re_var = 0)
      alpha = .22
      
      #get distribution of observed gammas
      glist = c()
      glist1 = c()
      glist2 = c()
      
      for (cl in unique(df$neigh)){
        #print('looping through clusters')
        subdf = df[df$neigh == cl,]
        if(length(unique(unlist(subdf[,trt_col]))) == 1){print('skip: no treatment variation'); next}
        if(length(unique(unlist(subdf[,cov_cols]))) == 1){print('skip: no covariate variation'); next}
        
        if(!bivar){ #univariate
          #print('estimating gamma range')
          mod = glm(unlist(subdf[,trt_col]) ~ (unlist(subdf[,cov_cols])), family = binomial()) #maybe need to add unlist if breaks
          glist = append(glist, coef(mod)[2])
        }
        
        if(bivar){ #bivariate
          mod1 = glm(unlist(subdf[,trt_col]) ~ (unlist(subdf[,cov_cols[1]])), family = binomial()) #maybe need to add unlist if breaks
          glist1 = append(glist1, coef(mod1)[2])
          
          mod2 = glm(unlist(subdf[,trt_col]) ~ (unlist(subdf[,cov_cols[2]])), family = binomial()) #maybe need to add unlist if breaks
          glist2 = append(glist2, coef(mod2)[2])
        }
        
      }
      
      if(!bivar){
        bounds = quantile(glist, probs = c(.1, .9))
        print('**********')
        print(i)
        print(ww)
        print(bounds)
        print('********')
        gamma_list = rbind(0,seq(from = bounds[1], 
                                 to = bounds[2], 
                                 length.out = gamma_R))
        gamma_list = cbind(gamma_list, c(0,0))
      }
      
      if(bivar){
        bounds1 = quantile(glist1, probs = c(.1, .9))
        bounds2 = quantile(glist2, probs = c(.1, .9))
        gl1 = c(seq(from = bounds1[1], to = bounds1[2], length.out = 7), 0)
        gl2 = c(seq(from = bounds2[1], to = bounds2[2], length.out = 7), 0)
        ggrid = expand.grid(gl1, gl2)
        gamma_list = t(cbind(0, ggrid))
      }
      
      gamma_numer = gamma_list
      ngam = ncol(gamma_list)
      
      yhat_group_sd2 = GroupIPW_sd2(dta = df, cov_cols = cov_cols,
                                    alpha = alpha, trt_col = trt_col, out_col = out_col,
                                    alpha_re_bound = 20, verbose = F,
                                    gamma_numer = gamma_list, loud_denom = F, fix_phi = T)
      
      z <- yhat_group_sd2$yhat_group
      po = apply(z, c(2,3), mean)
      
      haj = make_hajek(yhat_group_sd2, df)
      
      if(FALSE){
        boots_est <- BootVar_sd(dta = df, B = 100, alpha = alpha, gamma_numer = gamma_list,
                                cov_cols = cov_cols, trt_col = trt_col, out_col = out_col,
                                phi_hat_true = NULL, verbose = TRUE, return_everything = TRUE)
      }
      
      ypop = Ypop_sd(ygroup = yhat_group_sd2, horvitzthompson = F, dta = df)
      #ypop bootvar: bootvariance for each ypop by gamma
      if(FALSE){ypop_bootvar = apply(boots_est$boots, c(1,2), var)}
      
      
      yplot = data.frame(gammas = gamma_list[2,],
                         y0 = ypop$ypop[1,],
                         y1 = ypop$ypop[2,]) %>%
        pivot_longer(y0:y1)
      
      if(!bivar){
        ggplot(data = yplot, aes(x = gammas, y = value, group = name, color = name)) + 
          geom_line() +
          labs(y = 'potential outcome', title = paste0('Potential Outcomes (haj) in Cai Data\n', target_name))# + 
        #ylim(0, 1)
        #ggsave(here('figures', fig_loc, paste0('haj_yhat_cai_', target_name, '.png')), height = 5, width = 5)
        
        direct = DE_sd(ypop = ypop$ypop, ypop_var = ypop$ypop_var, boots = NULL)#boots = boots_est$hajboots)
        de_test = oe_sigtest(direct, ypop$ypop_var_de, seed = seed)
        
        de_univar_pvals = rbind(de_univar_pvals, 
                                c(i, 
                                  names(df)[i], 
                                  counter, 
                                  ww, 
                                  de_test$accept2, 
                                  de_test$pvalue))
        
        detab = data.frame(estimated.de = direct['est',],
                           lb = direct['est',] - 1.96*sqrt(direct['var',]),
                           ub = direct['est',] + 1.96*sqrt(direct['var',]),
                           gamma = gamma_list[2,],
                           target = target_name) 
        de_results[[counter]] <- detab
        ggplot(detab,
               aes(x = gamma, y = estimated.de)) + 
          geom_line() + 
          geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), alpha = .3) + 
          ggtitle(paste0('Direct Effect in Cai data \nReject null:',!de_test$accept, ' Pvalue = ', de_test$pvalue)) + 
          ylab('direct effect') + 
          labs(caption = 'range of gammas is the 10th and 90th percentile of cluster level trt~degree betas')
        if(revisions == T){
          ggsave(here('figures',fig_loc, paste0('haj_de_cai_', target_name, '_gammacount', gamma_R, '.png')), height = 5, width = 5)
        }else{
          ggsave(here('figures',fig_loc, paste0('haj_de_cai_', target_name, '.png')), height = 5, width = 5)
        }
        
        
        #indirect effect
        indirect0 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,1,], 
                          hajofygroup = haj$haj[1,], boots = NULL,#boots_est$hajboots, 
                          treatment = 1)
        ie0_cov = indirect0$ie_var
        indirect0 = indirect0$ie[,,ngam]
        ie0_test = oe_sigtest(indirect0, ie0_cov, seed = seed)
        
        ie0_univar_pvals = rbind(ie0_univar_pvals, 
                                 c(i, 
                                   names(df)[i], 
                                   counter, 
                                   ww, 
                                   ie0_test$accept2, 
                                   ie0_test$pvalue))
        
        
        indirect1 = IE_sd(ypop = ypop, ygroup = yhat_group_sd2$yhat_group[,2,], 
                          hajofygroup = haj$haj[2,], boots = NULL,#boots_est$hajboots, 
                          treatment =2)
        ie1_cov = indirect1$ie_var
        indirect1 = indirect1$ie[,,ngam]
        ie1_test = oe_sigtest(indirect1, ie1_cov, seed = seed)
        
        ie1_univar_pvals = rbind(ie1_univar_pvals, 
                                 c(i, 
                                   names(df)[i], 
                                   counter, 
                                   ww, 
                                   ie1_test$accept2, 
                                   ie1_test$pvalue))
        
        ievartab1 = data.frame(estimated.ie = indirect1['est',],
                               lb = indirect1['est',] - 1.96*sqrt(indirect1['var',]),
                               ub = indirect1['est',] + 1.96*sqrt(indirect1['var',]),
                               gamma = gamma_list[2,],
                               target = target_name) 
        ie1_results[[counter]] <- ievartab1
        
        ggplot(ievartab1,
               aes(x = gamma, y = estimated.ie)) + 
          geom_line() + 
          geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), alpha = .3) + 
          ggtitle(paste0('Indirect Effect in Cai data (on treated) \nReject null:',!ie1_test$accept, ' Pvalue = ', ie1_test$pvalue)) + 
          ylab('indirect effect') + 
          labs(caption = 'range of gammas is the 10th and 90th percentile of cluster level trt~degree betas')
        if(revisions == T){
          ggsave(here('figures',fig_loc, paste0('haj_ie1_cai_', target_name, '_gammacount', gamma_R, '.png')), height = 5, width = 5)
          
        }else{
          ggsave(here('figures',fig_loc, paste0('haj_ie1_cai_', target_name, '.png')), height = 5, width = 5)
        }
        
        ievartab0 = data.frame(estimated.ie = indirect0['est',],
                               lb = indirect0['est',] - 1.96*sqrt(indirect0['var',]),
                               ub = indirect0['est',] + 1.96*sqrt(indirect0['var',]),
                               gamma = gamma_list[2,],
                               target = target_name) 
        ie0_results[[counter]] <- ievartab0
        
        
        ggplot(ievartab0,
               aes(x = gamma, y = estimated.ie)) + 
          geom_line() + 
          geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub), alpha = .3) + 
          ggtitle(paste0('Indirect Effect in Cai data (on untreated) \nReject null:',!ie0_test$accept, ' Pvalue = ', ie0_test$pvalue)) + 
          ylab('indirect effect') + 
          labs(caption = 'range of gammas is the 10th and 90th percentile of cluster level trt~degree betas')
        #ggsave(here('figures',fig_loc, paste0('haj_ie1_cai_', target_name, '.png')), height = 5, width = 5)
        
        
        #overall effect
        oe = OE_sd(ypop = ypop, ygroup = yhat_group_sd2$oe_yhat_group,  
                   hajofygroup = haj$oe_haj, boots = NULL)#boots_est$oehajboots)
        oe_cov = oe$oe_var
        oe = oe$oe[,,ngam]
        oe_test = oe_sigtest(oe, oe_cov, seed = seed)
        
        print(list(i, 
                   names(df)[i], 
                   counter, 
                   ww, 
                   oe_test$accept2, 
                   oe_test$pvalue))
        univar_pvals = rbind(univar_pvals, c(i, 
                                             names(df)[i], 
                                             counter, 
                                             ww, 
                                             oe_test$accept2, 
                                             oe_test$pvalue))
        
        oevartab = data.frame(estimated.oe = oe['est',],
                              lb = oe['est',] - 1.96*sqrt(oe['var',]),
                              ub = oe['est',] + 1.96*sqrt(oe['var',]),
                              gamma = gamma_list[2,],
                              target = target_name) 
        oe_results[[counter]] <- oevartab
        
        ggplot(oevartab,
               aes(x = gamma, y = estimated.oe)) + 
          geom_line() + 
          theme_minimal() +
          theme(strip.text = element_text(size = 12),
                text = element_text(size = 12),
                axis.text = element_text(size = 12),
                legend.title = element_blank(),
                legend.position = 'bottom',
                panel.spacing = unit(1.4, "lines")) +
          geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub),alpha = 0.2,
                      fill = '#00BFC4', colour = rgb(0,0,0,0)) +
          ggtitle(paste0('Intervening on ', clean_name)) + 
          ylab('Overall Effect') + xlab(TeX(r'(\gamma)')) +
          labs(caption = paste0('Statistically significant deviation from flat OE? ',!oe_test$accept, ' Pvalue = ', oe_test$pvalue))
        if(revisions == T){
          ggsave(here('figures',fig_loc, paste0('haj_oe_cai_', paste0(target_name), '_gammacount', gamma_R,'.png')), height = 3, width = 5)
        }else{
          ggsave(here('figures',fig_loc, paste0('haj_oe_cai_', paste0(target_name), '.png')), height = 3, width = 5)
        }
        newline = c(target_name, cov_cols, gamma_R, oe_test$pvalue)
        if(!exists("oe_test_gamma_range_results")){
          oe_test_gamma_range_results = data.frame(
            'target_name' = c(),
            'covar' = c(), 
            'number_of_gammas' = c(),
            'pvalue' = c())}
        oe_test_gamma_range_results = rbind(oe_test_gamma_range_results, newline)
      }
      
      if(bivar){
        oe = OE_sd(ypop = ypop, ygroup = yhat_group_sd2$oe_yhat_group,  
                   hajofygroup = haj$oe_haj, boots = NULL)#boots_est$oehajboots)
        oe_cov = oe$oe_var
        oe = oe$oe[,,ngam]
        oe_test = oe_sigtest(oe, oe_cov, seed = seed)
        
        
        bivar_pvals = rbind(bivar_pvals, c(i, 
                                           names(df)[i], 
                                           counter, 
                                           ww, 
                                           oe_test$accept2, 
                                           oe_test$pvalue))
        
        #get list of scenarios where gamma is min for degree.
        #then only include those scenarios when testing for heterogeneity in effect
        use_cols = gamma_list[2,]==min(gamma_list[2,])
        oe_test_univar_bivar = oe_sigtest(oe[,use_cols], oe_cov[use_cols,use_cols], gammas_idx = c(1:sum(use_cols)), seed = seed)
        
        sink(file = paste0('test_bivar_univar/mindgr', ww ,'_', cov_cols[2],".txt"))
        print(cov_cols)
        print(ww)
        print(oe_test_univar_bivar$accept2)
        print(oe_test_univar_bivar$pvalue)
        sink(file = NULL)
        
        oevartab = data.frame(estimated.oe = oe['est',],
                              lb = oe['est',] - 1.96*sqrt(oe['var',]),
                              ub = oe['est',] + 1.96*sqrt(oe['var',]),
                              g1 = gamma_list[2,], #need to round these or make it so the label only shows 2 digits
                              g2 = gamma_list[3,],
                              target = target_name) 
        bivar_oe_results[[counter]] <- oevartab
        ggplot(oevartab,
               aes(x = g1, y = g2, fill = estimated.oe)) + 
          geom_raster() + 
          scale_fill_viridis() + 
          xlab(expression(paste(gamma['degree']))) +
          ylab(expression(paste(gamma['rice area']))) + 
          labs(fill = 'OE',
               #title = 'Overall Effect In Cai Data',
               caption = paste0('Statistically Significant? ',!oe_test$accept, ' Pvalue = ', oe_test$pvalue))
        if(revisions == T){
          ggsave(here('figures',fig_loc, paste0('haj_oe_cai_', 'bivardgrrice_', target_name, '_gammacount', gamma_R, '.png')), height = 5, width = 5)
        }else{
          ggsave(here('figures',fig_loc, paste0('haj_oe_cai_', 'bivardgrrice_', target_name, '.png')), height = 5, width = 5)
        }
        
        newline = data.frame(
          target_name = target_name,
          covar = j,
          number_of_gammas = gamma_R,
          pvalue = oe_test$pvalue
        )
        if(!exists("oe_test_gamma_range_results")){
          oe_test_gamma_range_results = data.frame(
            'target_name' = c(),
            'covar' = c(), 
            'number_of_gammas' = c(),
            'pvalue' = c())}
        oe_test_gamma_range_results = rbind(oe_test_gamma_range_results, newline)
        
      }
    }
    
  }
}

first_row <- gsub("^X\\.|\\.$", "", colnames(oe_test_gamma_range_results))

oe_test_gamma_range_results <- rbind(first_row, oe_test_gamma_range_results)

# update column names
colnames(oe_test_gamma_range_results) <- c("variable", "variable_num", "R_k", "pvalue")

write.csv(oe_test_gamma_range_results %>% arrange(variable), 
          file = '~/project_pi_lf474/sgd37/interference/figures/revision_figs_application/revisions_pval.csv')

if(bivar == F){
  print(univar_pvals)
  univar_pvals = data.frame(univar_pvals)
  namess = c('col_num', 'var_name', 'size_num', 'size_name', 'sig', 'pval')
  names(univar_pvals) = namess
  univar_pvals %>%  select(var_name, size_name, pval) %>% pivot_wider(names_from = size_name,
                                                                      values_from = c(pval)) %>% select(var_name, `_all`, `_le80`, `_over80`)
  
  de_univar_pvals = data.frame(de_univar_pvals)
  names(de_univar_pvals) = namess
  x = de_univar_pvals %>%  
    select(var_name, size_name, pval) %>% 
    pivot_wider(names_from = size_name,
                values_from = c(pval)) %>% 
    select(var_name, `_all`, `_le80`, `_over80`)
  
  ie0_univar_pvals = data.frame(ie0_univar_pvals)
  names(ie0_univar_pvals) = namess
  ie0_univar_pvals %>%  
    select(var_name, size_name, pval) %>% 
    pivot_wider(names_from = size_name,
                values_from = c(pval)) %>% 
    select(var_name, `_all`, `_le80`, `_over80`)
  
  ie1_univar_pvals = data.frame(ie1_univar_pvals)
  names(ie1_univar_pvals) = namess
  x = ie1_univar_pvals %>%  
    select(var_name, size_name, pval) %>% 
    pivot_wider(names_from = size_name,
                values_from = c(pval)) %>% 
    select(var_name, `_all`, `_le80`, `_over80`)
  
  
}
if(bivar){
  bivar_pvals = data.frame(bivar_pvals)
  names(bivar_pvals) = c('colnum1', 'colnum2', 'varname1', 'varname2', 'sizenum', 'sizename', 'sig', 'pval')
  bivar_pvals %>% select(varname2, sizename,pval) %>% 
    pivot_wider(names_from = sizename, values_from = pval) %>%
    select(varname2, `_all`, `_le80`, `_over80`)
}
######################################################
# COMP PLOTS #########################################
######################################################
which = c('Direct Effect', 'Indirect Effect (Treated)', 'Indirect Effect (Untreated)', 'Overall Effect')
counter = 0
for(res in list(de_results, ie1_results, ie0_results, oe_results)){ #
  #res = oe_results #counter = 4
  counter = counter+1
  estim_name = which[counter]
  de_results_all = list.rbind(res)
  #de_results_all[c('target1', 'target2', 'target3')] <- str_split_fixed(de_results_all$target, '_', n = 3)
  #de_results_all$target_var = ifelse(de_results_all$target3 == '', de_results_all$target1, paste0(de_results_all$target1, de_results_all$target2))
  #de_results_all$target_pop = ifelse(de_results_all$target3 == '', de_results_all$target2, de_results_all$target3)
  
  de_results_all$target_pop = sapply(str_split(de_results_all$target, '_'), function(x) tail(x,1))
  de_results_all$target_var = sapply(str_split(de_results_all$target, '_'), function(x) paste(head(x, -1), collapse = "_"))
  
  
  
  names(de_results_all)[1] = 'estimator'
  de_results_all = de_results_all %>% 
    mutate(target_pop = case_when(target_pop == 'all' ~ 'All Clusters',
                                  target_pop == 'le80' ~ '<= 80 units',
                                  target_pop == 'over80' ~ '>80 units')) %>%
    mutate(target_pop = factor(target_pop, 
                               levels = c('All Clusters', '<= 80 units' ,'>80 units')))
  
  
  
  for (i_var in unique(de_results_all$target_var)){
    #i_var = 'ricearea2010'
    sub = de_results_all %>% filter(target_var == i_var)
    
    clean_name2 = case_when(i_var == 'btwn' ~ 'Betweenness',
                            i_var == 'dgr' ~ 'Degree',
                            i_var == 'disaster_prob' ~ 'Future Disaster Probability',
                            i_var == 'ricearea_2010_imp' ~ 'Rice Production Area')
    
    ggplot(sub,
           aes(x = gamma, y = estimator)) + 
      geom_line() + 
      theme_minimal() +
      theme(strip.text = element_text(size = 12),
            text = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            legend.title = element_blank(),
            legend.position = 'bottom',
            panel.spacing = unit(1.4, "lines")) +
      geom_ribbon(aes(x = gamma, ymin = lb, ymax = ub),alpha = 0.2,
                  fill = '#00BFC4', colour = rgb(0,0,0,0)) +
      facet_wrap(~target_pop, scales = 'free_x') + 
      #ggtitle(paste0(estim_name, ' in Cai data')) + 
      ylab(estim_name) +  xlab(TeX(r'(\gamma)')) #+ 
    #ggtitle(paste0('Intervening on ', clean_name2))
    
    #labs(caption = 'range of gammas is the 10th and 90th percentile of cluster level trt~degree betas') +
    
    ggsave(here('figures',fig_loc, paste0(estim_name, i_var, '.png')), height = 3, width = 7)
    
  }
  
}



theme_set(theme_minimal()+
            theme(strip.text = element_text(size = 10),
                  text = element_text(size = 11),
                  axis.text = element_text(size = 11),
                  legend.title = element_blank(),
                  legend.position = 'bottom',
                  panel.spacing = unit(.5, "lines"),
                  strip.background = element_rect(color = rgb(0,0,0,0))))
#bivariate case
bivar_oe_results_all = list.rbind(bivar_oe_results) #%>%
#mutate(gamma_dgr = as.numeric(as.character(gamma_dgr)),
#       gamma_rice = as.numeric(as.character(gamma_rice)))

#clean labels for the clusters targeted
bivar_oe_results_all$target_pop = sapply(str_split(bivar_oe_results_all$target, '_'), function(x) tail(x,1))
bivar_oe_results_all$target_var = sapply(str_split(bivar_oe_results_all$target, '_'), function(x) paste(head(x, -1), collapse = "_"))

bivar_oe_results_all = bivar_oe_results_all %>% 
  mutate(target_pop = case_when(target_pop == 'all' ~ 'All Clusters',
                                target_pop == 'le80' ~ '<= 80 units',
                                target_pop == 'over80' ~ '>80 units'),
         target_pop = factor(target_pop, levels = c('All Clusters','<= 80 units','>80 units' ))) %>%
  filter(g1 !=0 & g2!=0) 

ggplot(bivar_oe_results_all,
       aes(x = g1, y = g2, fill = estimated.oe)) + 
  geom_raster() + 
  scale_fill_viridis() + 
  xlab(expression(paste(gamma['degree']))) +
  ylab(expression(paste(gamma['rice area']))) + 
  labs(fill = 'OE') + 
  facet_wrap(~target_pop,scales = 'free') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text()) 
#ggtitle(paste0(estim_name, ' in Cai data'))  +  xlab(TeX(r'(\gamma)'))
ggsave(here('figures', fig_loc, 'degreerice.png'), height = 3.5, width = 7)


test = bivar_oe_results[[3]]
ggplot(test, aes(gamma_dgr, gamma_rice)) + geom_point()







######################################################
# CORR PLOTS #########################################
######################################################


#correlation map of degree with other things
library(corrplot)
covs = c('educ', 'n.peers', 'age', 'rice_inc', 'ricearea_2010', 'disaster_loss', 'insurance_repay', 'insurance_buy',
         'disaster_prob', 'general_trust', 'default', 'info_takeup_rate', 'risk_averse', 'disaster_yes', 
         'literacy', 'educ_good', 'day', 'understanding',
         'cai_degree')

M= cor(df[,c(covs, 'outcome')], use ='pairwise.complete.obs')
corrplot(M)
ggsave(here('figures',fig_loc, paste0('corrplot', '.png')))


#among the treated, how much does degree effect the outcome?
treated = df %>% filter(trt == 1)
tmod = glm(data = treated, outcome ~ cai_degree, family = binomial())
summary(tmod) # not signif

tmod = glm(data = treated, outcome ~ cai_degree + cai_degree^2, family = binomial())
summary(tmod)

library(splines)
tmod <- glm(outcome ~ ns(cai_degree,4), data = treated, family = "poisson")
summary(tmod)

#how does the probability of outcome after treatment vary by degree
degreetab = df %>% group_by(outcome, cai_degree) %>% summarise(count=n()) %>%
  merge(df %>% group_by(cai_degree) %>% summarise(totalcount = n(), by = 'cai_degree')) %>%
  mutate(pct_uptake = count / totalcount, outcome = factor(outcome)) 

ggplot(data = degreetab %>% filter(outcome ==1), aes(x = cai_degree, y = pct_uptake)) + 
  geom_bar(stat='identity')
ggsave('feb16meeting/percentadoptbydegree.pdf')

ggplot(data = degreetab, aes(x = cai_degree, y = count, fill = outcome)) + 
  geom_bar(stat = 'identity')
ggsave('feb16meeting/countadoptbydegree.pdf')

#repeat for treated
degreetab = df %>% filter(trt ==1) %>% group_by(outcome, cai_degree) %>% summarise(count=n()) %>%
  merge(df %>% group_by(cai_degree) %>% summarise(totalcount = n(), by = 'cai_degree')) %>%
  mutate(pct_uptake = count / totalcount, outcome = factor(outcome)) 

ggplot(data = degreetab %>% filter(outcome ==1), aes(x = cai_degree, y = pct_uptake)) + 
  geom_bar(stat='identity')
ggsave('feb16meeting/treatedpercentadoptbydegree.pdf')

ggplot(data = degreetab, aes(x = cai_degree, y = count, fill = outcome)) + 
  geom_bar(stat = 'identity')
ggsave('feb16meeting/treatedcountadoptbydegree.pdf')

#exploratory: regress avg outcome of neighbors on other vars
library(igraph)
g = graph_from_adjacency_matrix(A)
outlist = c()
for (ind in V(g)){
  ne = neighbors(g, v = ind, mode = 'all')
  if(length(ne) > 0){
    den = length(ne)
    num = sum(df$outcome[df$id %in% names(ne)])
    outlist = append(outlist, num / den)
  }else{
    outlist = append(outlist, NA)
  }
}
outdf = data.frame(id = names(V(g)), out = outlist) %>% merge(df, by = 'id')
summary(outdf$out)

hist(outdf$out)

covs = names(outdf)[-c(1,2,3,4,5,6,9,10,13,24, 25,33,36, 37, 38)]
ff = as.formula(paste0('out ~ ', paste(covs, collapse = '+')))
mod = lm(data = outdf %>% filter(trt == 1), ff)
summary(mod)

co = cor(outdf[,c('out', covs)], use = 'pairwise.complete.obs')
corrplot(co)
