# Effective treatment allocation strategies under partial interference
This code accompanies the manuscript https://arxiv.org/abs/2504.07305, by Samantha G Dean, Georgia Papadogeorgou, and Laura Forastiere.

# Overview
Users can run model_sim.R and diffusion.sim to recreate the simulated results from the manuscript, or cai_application.R to recreate the application results from the manuscript. 

# Script Directory
CalcHeterogeneousTrueIE.R
- creates get_het_ie(), called in model_sim. 
- calculates the true value of the estimands that the proposed estimators are evaluated compared to

compile_diffusion.R
- loads in the simulation results from the diffusion scenario, calculated the summary statistics, and creates the output visualizations

compile_helper_funcs.R
- functions that help with summarising the simulation results, called within the compile_diffusion or compileACIC scripts

compileACIC.R

de_sd.R
- creates DE_sd(), called in model_sim
- calculated the direct effect from the input Y's

denominator_sd.R
- calculates the denominator of the proposed estiamators.

diffusion_sim.R

diffusion_truth_sim.R

diffusion_truth.R

GroupIPW_sd v2.R
- creates GroupIPW_sd2(), called in model_sim.R
- calculates the average outcome by cluster for the hypothetical treatment scenarios

hajek_adj.R

helper_functs.R

ie_sd.R

load_clean_cai.R

oe_sd.R

oe_stattest.R

parallel_bootvar_function_Test.R
- creates BootVar_sd(), called in model_sim.R

ypop_sd.R