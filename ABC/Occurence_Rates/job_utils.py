#This file is used to calculate the bin edges based on a number of bins/jobs

import numpy as np
import os

#Gets job_id from environment
job_id = int(os.getenv('JOB_ID'))

#Assign number of p and r bins
nbins = 100
num_rbins_per_run = 10
#Create arrays 
pmin, pmax = 3, 30
rmin, rmax = 1.4, 2.4

pmin_log, pmax_log = np.log(pmin), np.log(pmax)
rmin_log, rmax_log = np.log(rmin), np.log(rmax)

pbins = np.linspace(pmin_log, pmax_log, nbins) #Make it one smaller to differentiate 
rbins = np.linspace(rmin_log, rmax_log, nbins + 1)

def calculate_bins(job_id, pbins, rbins, nbins, num_rbins_per_run):
    '''
    This function calculates the period and rbin edges based on the number of total bins, rbins per run and job_id. 

    Input:
        job_id: int, variable set by user
        nbins: int, number of total bins
        num_rbins_per_run: int, number of rbins per run

    Output:
        pbins: np.array, period bin edges
        rbins: np.array, radius bin edges
        num_rbins_per_run: int, number of rbins per run
        pbin_lower: float, lower period bin edge of current job_id
        pbin_upper: float, upper period bin edge of current job_id
        r_start: int, start of rbins
        r_end: int, end of rbins
    '''

    
    #Compute the period and radius bin indices
    pbin_index = (job_id - 1) // (nbins // num_rbins_per_run) #Each 10 bins is 1 period bin
    rbin_index = (job_id - 1) % (nbins // num_rbins_per_run) #Cycles every 10 bins 
    
    #Define the period and radius bin edges
    pbin_lower, pbin_upper = pbins[pbin_index], pbins[pbin_index+1]
    r_start = rbin_index * num_rbins_per_run
    r_end = r_start + num_rbins_per_run
    
    return pbin_index, rbin_index, pbin_lower, pbin_upper, r_start, r_end

def job_info(nbins, num_rbins_per_run):
    '''
    This function creates a dictionary with the job_id, nbins, num_rbins_per_run, rbins, pbins, pbin_lower, 
    pbin_upper, r_start, r_end

    Input:
        nbins: int, number of total bins
        num_rbins_per_run: int, number of rbins per run

    Output:
        job_data: dict, dictionary containing the job_id, nbins, num_rbins_per_run, rbins, pbins, pbin_lower, 
        pbin_upper, r_start, r_end
    '''

    #Calculates bins
    pbin_index, rbin_index, pbin_lower, pbin_upper, r_start, r_end = calculate_bins(job_id, pbins, rbins, nbins, 
                                                                                             num_rbins_per_run)
    
    #Create dictionary 
    job_data = {
    "job_id": job_id,
    "nbins": nbins,
    "num_rbins_per_run": num_rbins_per_run,
    "rbins": rbins,
    "pbins": pbins,
    "pbin_index": pbin_index,
    "rbin_index": rbin_index,
    "pbin_lower": pbin_lower,
    "pbin_upper": pbin_upper,
    "r_start": r_start,
    "r_end": r_end,
    }
    return job_data