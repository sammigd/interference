# run this to compile results after running sbatch power_sims.sh
library(tidyverse)
library(xtable)

fig_loc <- "/nfs/roberts/project/pi_lf474/sgd37/interference/figures/ms_figs_univar_stdalpha_round2"

#load in files
oe_files    <- list.files(fig_loc, pattern = "^power_pset_\\d+\\.RData$",    full.names = TRUE)
ie_de_files <- list.files(fig_loc, pattern = "^power_ie_de_pset_\\d+\\.RData$", full.names = TRUE)

combine_rdata <- function(files, obj_name) {
  do.call(rbind, lapply(files, function(f) {
    e <- new.env()
    load(f, envir = e)
    get(obj_name, envir = e)
  }))
}

empty       <- combine_rdata(oe_files,    "empty")
empty_other <- combine_rdata(ie_de_files, "empty_other")

save(empty_other, file = file.path(fig_loc, 'bivar_power_ie_de.RData'))
save(empty,       file = file.path(fig_loc, 'bivar_power.RData'))

#create summary table
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
         rej_rate = `FALSE` / total,
         parms_copy = Var1) %>%
  separate(Var1, into = c('x', 'b3', 'b4', 'concordance', 'b5', 'ngam'), sep = '_') %>%
  mutate(ngam = str_sub(ngam, end = -7))
save(summary_tab, file = paste0(fig_loc, '/univar_power.RData') )

power_tab = teststattable %>% 
  merge(summary_tab, by.x = 'parms', by.y = 'parms_copy')  %>%
  dplyr::select(-x, -`TRUE`, -`FALSE`, -total) %>%
  arrange(b4, b5, b3, concordance, as.numeric(ngam)) %>%
  mutate(Level = ifelse(b4==0 & b5==0, round(rej_rate,2), ' '),
         Power = ifelse(b4!=0 | b5 !=0, round(rej_rate,2), ' ')) %>%
  dplyr::select(b3, b4, b5, concordance, ngam, teststat, Level, Power)

format_range <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return('')
  mn <- round(min(x), 2); mx <- round(max(x), 2)
  paste0(mn, ", ", mx)
}

power_tab_collapsed <- power_tab %>%
  group_by(b3, b4, b5, concordance) %>%
  summarise(
    #ngam     = format_range(ngam),
    teststat = format_range(teststat),
    Level    = format_range(Level),
    Power    = format_range(Power),
    .groups  = "drop"
  )

print(xtable(power_tab_collapsed, type = "latex"), file = paste0(fig_loc, "/ngam_power_table.tex"), include.rownames = F)

#end overall effect

#start ie and de

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
  separate(parms, into = c('x', 'b3', 'b4', 'concordance', 'b5', 'ngam'), sep = '_') %>%
  mutate(ngam = str_sub(ngam, end = -7),
         ngam = round(as.numeric(ngam)), 0) %>%
  arrange(b4, b5, b3, concordance, ngam) %>%
  rename(de_level = de_rej_rate) %>%
  mutate(Ie0_Level = ifelse(b4==0 & b5==0, round(ie0_rej_rate,2), ' '),
         Ie0_Power = ifelse(b4!=0 | b5 !=0, round(ie0_rej_rate,2), ' '),
         Ie1_Level = ifelse(b4==0 & b5==0, round(ie1_rej_rate,2), ' '),
         Ie1_Power = ifelse(b4!=0 | b5 !=0, round(ie1_rej_rate,2), ' '),
         de_level = round(de_level, 2), ngam = as.character(ngam)) %>%
  dplyr::select(b3, b4, b5, concordance, ngam, de_level, Ie0_Level, Ie0_Power, Ie1_Level, Ie1_Power)

other_tab_collapsed = other_power %>%
  group_by(b3, b4, b5, concordance) %>%
  summarise(
    #ngam     = format_range(ngam),
    de_level = format_range(de_level),
    Ie0_Level    = format_range(Ie0_Level),
    Ie0_Power    = format_range(Ie0_Power),
    Ie1_Level    = format_range(Ie1_Level),
    Ie1_Power    = format_range(Ie1_Power),
    .groups  = "drop"
  )

print(xtable(other_tab_collapsed, type = "latex"), file = paste0(fig_loc, '/power_table_ie_de.tex'), include.rownames = F)
