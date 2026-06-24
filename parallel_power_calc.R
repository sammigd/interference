#compile big sim results
library(tidyverse)
library(viridis)
library(ggh4x)
library(reporter)
library(xtable)
library(here)
library(ggplot2)
library(latex2exp)

###########################################################
# INPUTS TO SELECT SCENARIO ###############################
###########################################################

bivar = F #SET INPUT TO SELECT SCENARIO

if(bivar == F){ #SET INPUT TO SELECT SCENARIO (comment out all but one)
  revision_scenarios = 'round2' #second round of revisions, with varied # gammas
  #revision_scenarios = 'og' # univariate - original mod results - alpha = 0.5, clust size = 15
  #revision_scenarios = 'smallalpha' # univariate - revision addition - alpha = 0.4, clust size = 15
  #revision_scenarios = 'largealpha' #univariate - revision addition - alpha = 0.6, clust size = 15
  #revision_scenarios = 'altclus' #univariate - revision addition - alpha = 0.5, clust size = half 10, half 15
}

diffusion = F #SET INPUT TO SELECT SCENARIO

##################################################
# END SCENARIO INPUTS TO MANUALLY SET ############
##################################################

setwd("~/project_pi_lf474/sgd37")
source('interference/analysis_scripts/compile_helper_funcs.R')

if(bivar){fig_loc = "~/project/cai/figures/ms_figs_bivar_two"}

if(!bivar & revision_scenarios == 'og'){fig_loc = "~/project/cai/figures/ms_figs_univar_stdalpha"}
if(!bivar & revision_scenarios == 'smallalpha'){fig_loc = "~/project/cai/figures/ms_figs_univar_smallalpha"}
if(!bivar & revision_scenarios == 'largealpha'){fig_loc = "~/project/cai/figures/ms_figs_univar_bigalpha"}
if(!bivar & revision_scenarios == 'altclus'){fig_loc = "~/project/cai/figures/ms_figs_univar_altclus"}
if(!bivar & revision_scenarios == 'round2'){fig_loc = "~/project_pi_lf474/sgd37/interference/figures/ms_figs_univar_stdalpha_round2"}

#univar heterogeneity
if(!bivar){
  make_gamma_table = function(ngam, third){
    #ngam = 100
    #third = 33
    gl = seq(from = -1.3, to = 1.3, length.out = third)
    ngl = rep(0, third)
    gamma_list = rbind(rep(0, third*3),
                       c(gl, ngl, ngl),
                       c(ngl, gl, ngl),
                       c(ngl, ngl, gl))
    gamma_list = cbind(gamma_list, c(0,0,0,0))
    return(gamma_list)
  }
  
  #gamma_list = make_gamma_table(100, 33)
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

#set folders to read in raw simulation results
if(bivar){setwd("~/project/cai/ms_bivar_2024jun01")}

if(!bivar & revision_scenarios == 'og'){setwd("~/project/cai/ms_univar_2025jul9_stdalpha")}
if(!bivar & revision_scenarios == 'smallalpha'){setwd("~/project/cai/ms_univar_2025jul9_smallalpha")}
if(!bivar & revision_scenarios == 'largealpha'){setwd("~/project/cai/ms_univar_2025jul9_largealpha")}
if(!bivar & revision_scenarios == 'altclus'){setwd("~/project/cai/ms_univar_2025aug17_altclus")}
if(!bivar & revision_scenarios == 'round2'){setwd("~/project_pi_lf474/sgd37/interference/results/ms_univar_2026may16_stdalpha")}

alliters = list.files()

# CHECK PARAMETERS INCLUDED IN SIMUALTIONS IN OUTPUT
table(str_remove(alliters, word(alliters, sep = '_')))
parmlist = unique(str_remove(alliters, word(alliters, sep = '_')))#[-25]
print(parmlist) # CHECK THAT ALL THE .RSAVE FILE NAMES LOOK RIGHT

source('/nfs/roberts/project/pi_lf474/sgd37/interference/analysis_scripts/oe_stattest.R')

#get counts
table(sapply((str_split(alliters, '_', n = 2)), tail , 1))
sigresults = c()
sigresults_ie1 = c()
sigresults_ie0 = c()
sigresults_de = c()
empty = matrix(nrow = length(alliters), ncol = 3)
empty_other = matrix(nrow = length(alliters), ncol = 4)

ind = 1

power_parmlist = parmlist
#for each set of parameters

args        <- commandArgs(trailingOnly = TRUE)
pset_index  <- as.integer(args[1])
revision    <- args[2]   # round2 label

# INITIAL PATHS
base_dir <- "/nfs/roberts/project/pi_lf474/sgd37"
setwd(base_dir)
source('interference/analysis_scripts/compile_helper_funcs.R')
source('interference/analysis_scripts/oe_stattest.R')

if (revision_scenarios == 'round2') {
  results_dir <- file.path(base_dir, "interference/results/ms_univar_2026may16_stdalpha")
  fig_loc     <- file.path(base_dir, "interference/figures/ms_figs_univar_stdalpha_round2")
}

#BUILD PARAMETER LIST
alliters      <- list.files(results_dir)
parmlist      <- unique(str_remove(alliters, word(alliters, sep = '_')))
power_parmlist <- parmlist

#SELECT PARAMETERS FOR THIS TASK
pset <- power_parmlist[pset_index]
iters <- alliters[str_detect(alliters, fixed(pset))]

message("Task ", pset_index, " | pset: ", pset, " | n_iters: ", length(iters))

#INITIATE OUTPUT OBJECT
n         <- length(iters)
empty     <- matrix(NA_character_, nrow = n, ncol = 3)
empty_other <- matrix(NA_character_, nrow = n, ncol = 4)

ind <- 1
for (i in seq_along(iters)) {
  load(file.path(results_dir, iters[i]))
  if (is.na(test$oe[1, 1])) next
  
  testresult     <- oe_sigtest(test$oe,        test$oe_cov,            1:dim(test$oe_cov)[1])
  testresult_ie1 <- oe_sigtest(test$indirect1, test$indirect_cov[,,2], 1:dim(test$oe_cov)[1])
  testresult_ie0 <- oe_sigtest(test$indirect0, test$indirect_cov[,,1], 1:dim(test$oe_cov)[1])
  testresult_de  <- oe_sigtest(test$direct,    test$direct_cov)
  
  empty[ind, ]       <- c(pset, testresult$accept,    testresult$diff)
  empty_other[ind, ] <- c(pset, testresult_ie0$accept, testresult_ie1$accept, testresult_de$accept)
  ind <- ind + 1
}

#SAVE RESULTS
out_file      <- file.path(fig_loc, paste0("power_pset_", pset_index, ".RData"))
out_file_other <- file.path(fig_loc, paste0("power_ie_de_pset_", pset_index, ".RData"))
save(empty,       file = out_file)
save(empty_other, file = out_file_other)

message("Done. Saved to ", out_file)
