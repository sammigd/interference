# Effective treatment allocation strategies under partial interference
This code accompanies the manuscript https://arxiv.org/abs/2504.07305, by 
Samantha G Dean, Georgia Papadogeorgou, and Laura Forastiere.

# Instructions
To replicate the simulation study
- Download this github repo
- Run simulation_scripts/model_sim.R. Specify in the file inputs if you want to 
  run bivariate/univariate, and linear model/diffusion settings. The file is 
  configured to be sent to run in parallel on an hpc. See example shell scripts 
  in same folder. When running from a cluster, make sure your current directory
  contains the model_sim.R script when you run the shell script. You will need 
  to update the working directory in the script to run on your computer.  
  
To replicate the application results
- Download this github repo
- Run the script "cai_application.R". Customize inputs for variables of interest
  and bivariate versus univariate analysis.

To replicate power calculations
- After generating simulation study results, run parallel_power_calc.R (recommend using sbatch power_sim.sh) and then compile_parallel_power.R


Other scripts, within the 'analysis_scripts' folder, create functions that are 
called within the model_sim script. These scripts are adapted from our 
coauthor's prior work (https://github.com/gpapadog/Interference).
