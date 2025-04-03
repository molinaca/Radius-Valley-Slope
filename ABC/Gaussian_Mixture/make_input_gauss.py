# This file is used to create an input file for the package cosmoabc. The input file changes based on the results
# of the job_utils.py file. 
#Documentation for how to set up the input file can be found at https://pypi.org/project/cosmoabc/

import numpy as np
import os

#######################################

## Create input file
f = open(f'input_file_gapfit.txt','w+')
f.write('path_to_obs = ./Exo-Occurrence/Data/FGK_planets.dat # \n')

#Write paramteres r1-r{number_rbins}

f.write('param_to_fit = n1 mean1p mean1r cov1p cov1pr cov1r n2 mean2p mean2r cov2p cov2pr cov2r')
f.write(' # \n')

f.write('param_to_sim = n1 mean1p mean1r cov1p cov1pr cov1r n2 mean2p mean2r cov2p cov2pr cov2r')
f.write(' # \n')

f.write('prior_func =')
#Define prior function
n_param = 12
b = 0
while b < n_param:
	f.write(' flat_prior')
	b += 1
f.write(' # \n')

#Get paramter values for first gaussian 
f.write(f'n1_prior_par_name = pmin pmax # \n')
f.write(f'n1_prior_par_val = 0.15 0.28 # \n')
f.write(f'n1_lim = 0.10 0.30 # \n')
f.write(f'n1 = 0.2 # \n')

f.write(f'mean1p_prior_par_name = pmin pmax # \n')
f.write(f'mean1p_prior_par_val = {np.log10(3)} {np.log10(10)} # \n')
f.write(f'mean1p_lim = {np.log10(2)} {np.log10(15)} # \n')
f.write(f'mean1p = {np.log10(5)} # \n')

f.write(f'mean1r_prior_par_name = pmin pmax # \n')
f.write(f'mean1r_prior_par_val = {np.log10(1.2)} {np.log10(1.7)} # \n')
f.write(f'mean1r_lim = {np.log10(1.0)} {np.log10(2.0)} # \n')
f.write(f'mean1r = {np.log10(1.4)} # \n')

f.write(f'cov1p_prior_par_name = pmin pmax # \n')
f.write(f'cov1p_prior_par_val = 0.025 0.1 # \n')
f.write(f'cov1p_lim = 0.025 0.1 # \n')
f.write(f'cov1p = 0.075 # \n')

f.write(f'cov1pr_prior_par_name = pmin pmax # \n')
f.write(f'cov1pr_prior_par_val = -0.01 0.01 # \n')
f.write(f'cov1pr_lim = -0.01 0.01 # \n')
f.write(f'cov1pr = -0.01 # \n')

f.write(f'cov1r_prior_par_name = pmin pmax # \n')
f.write(f'cov1r_prior_par_val = 0.01 0.08 # \n')
f.write(f'cov1r_lim = 0.01 0.08 # \n')
f.write(f'cov1r = 0.02 # \n')

#Get paramtere values for second gaussian

f.write(f'n2_prior_par_name = pmin pmax # \n')
f.write(f'n2_prior_par_val = 0.15 0.28 # \n')
f.write(f'n2_lim = 0.10 0.30 # \n')
f.write(f'n2 = 0.2 # \n')

f.write(f'mean2p_prior_par_name = pmin pmax # \n')
f.write(f'mean2p_prior_par_val = {np.log10(8)} {np.log10(20)} # \n')
f.write(f'mean2p_lim = {np.log10(5)} {np.log10(30)} # \n')
f.write(f'mean2p = {np.log10(15)} # \n')

f.write(f'mean2r_prior_par_name = pmin pmax # \n')
f.write(f'mean2r_prior_par_val = {np.log10(2.2)} {np.log10(2.7)} # \n')
f.write(f'mean2r_lim = {np.log10(2.0)} {np.log10(2.8)} # \n')
f.write(f'mean2r = {np.log10(2.4)} # \n')

f.write(f'cov2p_prior_par_name = pmin pmax # \n')
f.write(f'cov2p_prior_par_val = 0.025 0.1 # \n')
f.write(f'cov2p_lim = 0.025 0.1 # \n')
f.write(f'cov2p = 0.075 # \n')

f.write(f'cov2pr_prior_par_name = pmin pmax # \n')
f.write(f'cov2pr_prior_par_val = -0.01 0.01 # \n')
f.write(f'cov2pr_lim = -0.01 0.01 # \n')
f.write(f'cov2pr = -0.01 # \n')

f.write(f'cov2r_prior_par_name = pmin pmax # \n')
f.write(f'cov2r_prior_par_val = 0.01 0.08 # \n')
f.write(f'cov2r_lim = 0.01 0.08 # \n')
f.write(f'cov2r = 0.02 # \n')


#Define cosmoabc parameters
f.write('M = 500 # \n') #Number of particles per system
f.write('Mini = 1000 # \n') #Number of possible Ndraws 
f.write('delta = 0.1 # \n') #Threshold for convergence
f.write('qthreshold = 0.75 # \n') #Quantile threshold
f.write(f'file_root = ./Gapfit_Output/rvalley_sims_1_ # \n') #Data file names
f.write('screen = 0 # \n') 
f.write('ncores = 1 # \n')
f.write('split_output = 1 # \n')
f.write('simulation_func = simulate_catalogue # \n') #Simulation function in function file
f.write('distance_func = distance # \n') #Distance function in function file
f.write('dist_dim = 1 # \n')
f.write('nruns = -9 # \n') #Number of runs, either defined or -9 to run until convergence
f.close()


