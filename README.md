# SCV_sim_yeast
Simulations and yeast data experiments for "On cross-validation for sparse reduced rank regression"

## Simulations
Simulated data experiments are in /simulations. Run paper_sim.m (main tables) and paper_sim_largep_largerb.m (extra table). Output will be Excel files with results. Edit the methods variable to turn on/off model aggregation.

## Applied data
Yeast cell cycle data experiments are in /yeast. Open yeast.Rproj (R project) then explore.R to create the data sets (install pacman package before running any of the R scripts). Experiments are run using the bootstrap_std_type2.m Matlab file - experiments can be broken up into chunks by editing the "curr_range" variable (line 132). Then results are processed using combine_bootstrap_std_type2.m. Plots are created using the bootstrap_std_type2_boxplots.R and bootstrap_std_type2_freq_plots.R files (open in the R project).
