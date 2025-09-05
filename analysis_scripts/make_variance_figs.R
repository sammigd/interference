library(tidyverse)
library(ggh4x)
library(latex2exp)

fig_loc = "~/project/cai/figures/"

load("~/project/cai/figures/ms_figs_univar_smallalpha/OE_figdf.RData")
smallalpha = fig2df
smallalpha$alpha = '0.4'
smallalpha$bias = NULL

load("~/project/cai/figures/ms_figs_univar_bigalpha/OE_figdf.RData")
bigalpha = fig2df
bigalpha$alpha = '0.6'

load("~/project/cai/figures/ms_figs_univar_stdalpha/OE_figdf.RData")
stdalpha = fig2df
stdalpha$alpha = '0.5'


varfigdf = rbind(smallalpha, bigalpha, stdalpha)

theme_set(theme_minimal()+
            theme(strip.text = element_text(size = 10),
                  text = element_text(size = 11),
                  axis.text = element_text(size = 11),
                  #legend.title = element_blank(),
                  legend.position = 'bottom',
                  panel.spacing = unit(.5, "lines"),
                  strip.background = element_rect(color = rgb(0,0,0,0))))

ggplot(varfigdf, aes(x = g, y = oe_bootvar_haj, colour = alpha)) + 
  geom_line() + 
  facet_nested(rows = vars(gamma_ind, Conc), cols = vars(IntIn1, Interference), nest_line = element_line(), solo_line = T) + 
  ylab('Bootstrapped Variance') + 
  xlab(TeX(r'(\gamma)')) + 
  labs(colour = 'hypothetical alpha')
ggsave(paste0(fig_loc, '/variancecompboot.png'), width = 9, height = 6)


ggplot(varfigdf, aes(x = g, y = oe_analyticalvar_haj, colour = alpha)) + 
  geom_line() + 
  facet_nested(rows = vars(gamma_ind, Conc), cols = vars(IntIn1, Interference), nest_line = element_line(), solo_line = T) + 
  ylab('Analytical Variance') + 
  xlab(TeX(r'(\gamma)')) + 
  labs(colour = 'hypothetical alpha')
ggsave(paste0(fig_loc, '/variancecompanalytical.png'), width = 9, height = 6)

