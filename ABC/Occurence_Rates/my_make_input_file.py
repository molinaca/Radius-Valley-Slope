# This file is used to create an input file for the package cosmoabc. The input file changes based on the results
# of the job_utils.py file. 
#Documentation for how to set up the input file can be found at https://pypi.org/project/cosmoabc/

import numpy as np
import os
from job_utils import job_info, nbins, num_rbins_per_run

#Use job_utils to gets bin information
job_data = job_info(nbins, num_rbins_per_run)

job_id = job_data['job_id']
pbins = job_data['pbins']
rbins = job_data['rbins']
pbin_index = job_data['pbin_index']
rbin_index = job_data['rbin_index']
pbin_lower = job_data['pbin_lower']
pbin_upper = job_data['pbin_upper']
r_start = job_data['r_start']

#######################################

## Log values of bins for further use
file = 'log.txt'

#If file does not exist create file
if not os.path.exists(file):
	with open(file, 'w') as f:
		f.write('')

#Read contents of files
with open(file, 'r') as f:
	file_contents = f.read()

#Create new line to append to file
p = np.exp(pbins[pbin_index])
new_line = f'Index: {job_id}, Pbin: {pbin_index}, {np.round(p,2)}'
for i in range(num_rbins_per_run):
	#Add rbin information to new line
	r = np.exp(rbins[r_start+i])
	new_line += f', Rbin{i+1}: {r_start+i}, {np.round(r, 2)}'
new_line += '\n'

#Write new info
with open(file, 'a') as f2:
	f2.write(new_line)

#######################################

## Create input file
f = open(f'input_file_{job_id}.txt','w+')
f.write('path_to_obs = /home/molinaca//Data/FGK_planets.dat # \n')

#Write paramteres r1-r{number_rbins}
b = 0
f.write('param_to_fit =')
while b < num_rbins_per_run:
	b += 1
	f.write(' r' + str(b))
f.write(' # \n')
b = 0
f.write('param_to_sim =')
while b < num_rbins_per_run:
	b += 1
	f.write(' r' + str(b))
f.write(' # \n')
b = 0
f.write('prior_func =')
#Define prior function
while b < num_rbins_per_run:
	b += 1
	f.write(' flat_prior')
f.write(' # \n')

#Calculate maximum occurence rate (fmax) for range of parameters 
b = 0
while b < num_rbins_per_run:
	f.write('r' + str(b+1) + '_prior_par_name = pmin pmax # \n')
	pmin, pmax = pbin_lower, pbin_upper
	rmin, rmax = rbins[r_start+b], rbins[r_start+(b+1)]
	fmax = 2*np.log2(np.exp(pmax-pmin))*np.log2(np.exp(rmax-rmin))
	f.write('r' + str(b+1) + '_prior_par_val = 0.0 ' + str(fmax) + ' # \n')
	f.write('r' + str(b+1) + '_lim = 0.0 ' + str(fmax) + ' # \n')
	f.write('r' + str(b+1) + ' = 0.05 # \n')
	b += 1

#Define cosmoabc parameters
f.write('M = 500 # \n') #Number of particles per system
f.write('Mini = 1000 # \n') #Number of possible Ndraws 
f.write('delta = 0.1# \n') #Threshold for convergence
f.write('qthreshold = 0.75 # \n') #Quantile threshold
f.write(f'file_root = ./Output/rvalley_sims_{job_id}_ # \n') #Data file names
f.write('screen = 0 # \n') 
f.write('ncores = 1 # \n')
f.write('split_output = 1 # \n')
f.write('simulation_func = simulate_catalogue # \n') #Simulation function in function file
f.write('distance_func = distance # \n') #Distance function in function file
f.write('dist_dim = 1 # \n')
f.write('nruns = -9 # \n') #Number of runs, either defined or -9 to run until convergence
f.close()


