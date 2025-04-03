#This file contains the functions used to calculate the distance between the simulated and observed data used 
# with cosmoabc.
#Functions taken from https://github.com/mkunimoto/Exo-Occurrence/blob/master/Codes/func_ExoOccurrence.py

import numpy as np
import sys
from scipy.stats import poisson
from scipy.stats import bernoulli
from scipy.stats import gamma
from scipy.stats import halfnorm
from job_utils import job_info, nbins, num_rbins_per_run

#Load star data
stars = np.loadtxt("./Exo-Occurrence/Data/FGK_properties.dat")
Nstars = len(stars)

##Constants##

#Define min and max 
rmin, rmax = 1.4, 2.4 #Rearth
pmin, pmax = 3, 30 #days

#Get bin information 
job_data = job_info(nbins, num_rbins_per_run)

job_id = job_data['job_id']
rbins = job_data['rbins']
pbin_lower = job_data['pbin_lower']
pbin_upper = job_data['pbin_upper']
r_start = job_data['r_start']

#Get data for transit calculations
earth_radius_in_m = 6.378e6
sun_radius_in_m = 6.9566e8
rsol_in_au = 0.0465
rearth_in_rsol = earth_radius_in_m/sun_radius_in_m

kepler_exp_time_internal = 6.019802903/(24.*60.*60.)
kepler_read_time_internal = 0.5189485261/(24.*60.*60.)
num_exposures_per_LC = 270.
LC_integration_time = kepler_exp_time_internal*num_exposures_per_LC
LC_read_time = kepler_read_time_internal*num_exposures_per_LC
LC_duration = LC_integration_time + LC_read_time
LC_rate = 1./LC_duration

##Functions##

def line(x, m, x0, y0):
	'''
	This function calculates the equation of the line.

	Input: 
		x: float, x value
		m: float, slope
		x0: float, x intercept
		y0: float, y intercept

	Output:
		y: float, y value
	'''
	return m*(x-x0)+y0

def remove_gap_planets(x, y, line, gap_width):

	'''
	This function removes planets that are in the radius valley gap. 

	Input:
		x: np.array, x values
		y: np.array, y values
		line: np.array, line values
		gap_width: float, width of the gap

	Output:
		x_new: np.array, x values without gap planets
		y_new: np.array, y values without gap planets
	'''
	
	dy = y - line
	remove = np.abs(dy) < gap_width
	x_new, y_new = [a[~remove] for a in (x, y)]
	
	return x_new, y_new

def semimajor_axis_radius_ratio(P, M, R):
	'''
	This function calculates the semimajor axis (a) to stellar radius (R) ratio.

	Input:
		P: float, period in days
		M: float, stellar mass in solar masses
		R: float, stellar radius in solar radii
	
	Output:
		a/R: float, semimajor axis to stellar radius ratio	
	'''
	a3 = M*(P/365.)**2. #Calculate a^3 
	a = a3**(1./3.)
	return a/R/rsol_in_au

def transit_depth(Rp, R, u1, u2):

	'''
	This function calculates the transit depth.

	Input:
		Rp: float, planet radius in Earth radii
		R: float, stellar radius in solar radii
		u1: float, limb darkening parameter
		u2: float, limb darkening parameter

	Output:
		depth: float, transit depth
	'''
	k = Rp/R*rearth_in_rsol
	c0 = 1. - (u1 + u2)
	w = c0/4. + (u1 + 2.*u2)/6. - u2/8.
	return 1. - (c0/4. + (u1+2.*u2)*(1.-k**2.)**(3./2.)/6. - u2*(1.-k**2)/8.)/w

def transit_duration(P, aR, b):
		
	'''
	This function calculates the transit duration.

	Input:
		P: float, period in days
		aR: float, semimajor axis to stellar radius ratio
		b: float, impact parameter
	
	Output:
		Tdur: float, transit duration
	'''

	Tdur = 24.*(P/np.pi)/aR*np.sqrt(1.-b**2.)
	
	return Tdur

def impact_parameter(aR, i):

	'''
	This function calculates the impact parameter.

	Input:
		aR: float, semimajor axis to stellar radius ratio
		i: float, inclination in radians
	
	Output:
		b: float, impact parameter
	'''

	b = aR*np.cos(i)

	return b

def num_transits(P, Tobs, fduty):

	'''
	This function calculates the number of transits.

	Input:
		P: float, period in days
		Tobs: float, observation time in days
		fduty: float, duty cycle

	Output:
		Ntr: float, number of transits
	'''

	Ntr = np.round(Tobs/P*fduty)

	return Ntr

def transit_SNR(P, Rp, dur, R, u1, u2, Tobs, fduty, std):

	'''
	This function calculates the transit signal-to-noise ratio.

	Input:
		P: float, period in days
		Rp: float, planet radius in Earth radii
		dur: float, transit duration in hours
		R: float, stellar radius in solar radii
		u1: float, limb darkening parameter
		u2: float, limb darkening parameter
		Tobs: float, observation time in days
		fduty: float, duty cycle
		std: float, standard deviation of noise

	Output:
		SNR: float, signal-to-noise ratio
	'''

	Ntr = num_transits(P, Tobs, fduty)
	depth = transit_depth(Rp, R, u1, u2)

	SNR = depth/std*np.sqrt(Ntr*dur*LC_rate)

	return SNR

def P_win(P, Tobs, fduty):

	'''
	This function calculates the window function.

	Input:
		P: float, period in days
		Tobs: float, observation time in days
		fduty: float, duty cycle

	Output:
		Pwin: float
	'''

	j = Tobs/P
	Pwin = np.where(j < 2., 0.0, 1.-(1.-fduty)**j-j*fduty*(1.-fduty)**(j-1.)-0.5*j*(j-1.)*fduty*fduty*(1.-fduty)**(j-2.))
	return Pwin

def det_comp(SNR, A, B, C):
		
	'''
	This function calculates the detection completeness.

	Input:
		SNR: float, signal-to-noise ratio
		A: float, parameter
		B: float, parameter
		C: float, parameter
	
	Output:
		det_comp: float, detection completeness
		
	'''
	det_comp = C*gamma.cdf(SNR, A, scale=B)
	
	return det_comp

def P_det(Ntr, SNR):

	'''
	This function calculates the detection probability.

	Input:
		Ntr: float, number of transits
		SNR: float, signal-to-noise ratio

	Output:
		Pdet: float, detection probability
	'''

	Pdet = np.where(Ntr >= 37, det_comp(SNR,12.23315232,  0.78203581,  0.94645662), np.where(Ntr >= 19, 
						det_comp(SNR,14.86511928,  0.72917663,  0.91642721), np.where(Ntr >= 10, 
						det_comp(SNR,11.45382365,  1.07249528,  0.91176524), np.where(Ntr >= 7, 
						det_comp(SNR,11.54128644,  1.20628098,  0.88169029), np.where(Ntr == 6, 
						det_comp(SNR,11.48116478,  1.27632116,  0.83694848), np.where(Ntr == 5, 
						det_comp(SNR,13.54878807,  1.09003313,  0.7704247), np.where(Ntr == 4, 
						det_comp(SNR,17.47440559,  0.90589395,  0.66508744), np.where(Ntr == 3, 
						det_comp(SNR,12.02912833,  1.38916308,  0.46525859), 0.0))))))))

	return Pdet

def transit_noise_model(P, Rp, R, aR, b, Tobs, fduty, std):

	'''
	This function calculates the transit noise model.

	Input:
		P: float, period in days
		Rp: float, planet radius in Earth radii
		R: float, stellar radius in solar radii
		aR: float, semimajor axis to stellar radius ratio
		b: float, impact parameter
		Tobs: float, observation time in days
		fduty: float, duty cycle
		std: float, standard deviation of noise

	Output:
		noise_model: float, transit noise model
	'''

	t0 = P*np.random.uniform()
	tao0 = P/(aR*2.*np.pi)
	r = rearth_in_rsol*Rp/R
	T = 2.*tao0*np.sqrt(1.-b**2.)
	tau = 2.*tao0*r/np.sqrt(1.-b**2.)
	Ttot = P
	I = LC_integration_time
	Lambda = LC_rate*num_transits(P, Tobs, fduty)
	tau3 = tau**3.
	I3 = I**3.
	a2 = (5.*tau3 + I3 - 5.*tau*tau*I)/tau3
	a3 = (9.*I**5.*Ttot - 40.*tau3*I*I*Ttot + 120.*tau**4.*I*(3.*Ttot - 2.*tau))/tau**6.
	a4 = (a3*tau**5. + I**4.*(54.*tau - 35.*Ttot) - 12.*tau*I3*(4.*tau+Ttot) + 360.*tau**4.*(tau - Ttot))/tau**5.
	a5 = (a2*(24.*T*T*(I-3.*tau) - 24.*T*Ttot*(I-3.*tau)) + tau3*a4)/tau3
	a11 = (I*Ttot - 3.*tau*(Ttot - 2.*tau))/(tau*tau)
	b1 = (6.*I*I - 3.*I*Ttot + tau*Ttot)/(I*I)
	b4 = (6.*T*T - 6.*T*Ttot + I*(5.*Ttot - 4.*I))/(I*I)
	b6 = (12.*b4*I3 + 4.*tau*(-6.*T*T + 6.*T*Ttot + I*(13.*Ttot - 30.*I)))/I3
	b7 = (b6*I**5. + 4.*tau*tau*I*I*(12.*I - 11.*Ttot) + tau3*I*(11.*Ttot - 6.*I) - tau**4.*Ttot)/I**5.
	depth = r**2.
	sigma_depth = np.where(tau >= I, np.sqrt(abs(24.*a11*a2/(tau*a5))), np.sqrt(abs(24.*b1/(I*b7))))
	sigma_depth *= std/np.sqrt(Lambda)
	noise_model = np.random.normal(loc=depth, scale=sigma_depth)
	return noise_model

def assumed_stellar_radius(sigdown, sigup, flag):
	'''
	This function calculates the assumed stellar radius.

	Input:
		sigdown: float, lower uncertainty in stellar radius
		sigup: float, upper uncertainty in stellar radius
		flag: int, flag for upper or lower uncertainty

	Output:
		R_assumed: float, assumed stellar radius
	'''

	R_assumed = np.where(flag == 1, halfnorm.rvs(scale=sigup), -halfnorm.rvs(scale=sigdown))

	return R_assumed

def simulate_catalogue(params):

	'''
	This function simulates a catalogue of planets and returns an array with the period and radius of each planet.

	Input:
		params: dict, dictionary of parameters

	Output:
		data_sim: np.array, array of simulated data
	'''

	#Get rates for each bin
	rates = [params[f'r{i}'] for i in range(1,num_rbins_per_run+1)]
	total_rate = np.sum(rates)
	prob_bin = rates/total_rate 
	num_bin = len(rates)

	#Simulate planet population 
	num_pl = poisson.rvs(total_rate*Nstars)

	#Assign radius bin to each planet
	if (num_bin == 1):
		Rp = np.exp(np.random.uniform(rbins[r_start], rbins[r_start+1],size=num_pl))
	else:
		bins = np.random.choice(num_bin, num_pl, p=prob_bin)
		Rp = np.exp(np.random.uniform(rbins[r_start+bins], rbins[r_start+(bins+1)]))
	P = np.exp(np.random.uniform(pbin_lower, pbin_upper,size=num_pl))
	# Get number of plaanets
	#total_rate = 0.16

	# Generate values for parameters
	#slope = np.exp(np.log(np.abs(params['slope'])))

	#if params['slope'] < 0:
	#	slope = -slope
    
	#x_offset = np.exp(np.log(params['x0']))
	#y_offset = np.exp(np.log(params['y0']))
	#gap_width = np.exp(np.log(params['width']))

	#xps = np.exp(np.random.uniform(np.log(xmin),np.log(xmax),size=num_pl))
	#Rps = np.exp(np.random.uniform(np.log(Rmin),np.log(Rmax),size=num_pl))

	# Create radius valley with parameters and remove planets in the gap
	#y = line(xps, slope, x_offset, y_offset)
	#xps_new, Rps_new = remove_gap_planets(xps, Rps, y, gap_width)

	#new_num_pl = len(xps_new)

    #Assign a random star to each planet

	#Assign a star to each planet
	allID = np.random.random_integers(0,high=len(stars)-1,size=num_pl)
	
	
    #Define stars columns:
	stars_chosen = stars[allID]
	sigdown = stars_chosen[:,2]
	sigup = stars_chosen[:,3]
	M_star = stars_chosen[:,5]
	u1 = stars_chosen[:,9]
	u2 = stars_chosen[:,10]
	std = stars_chosen[:,11]
	Tobs = stars_chosen[:,12]
	fduty = stars_chosen[:,13]

	# Calculate orbital properties based on stars 
	alli = np.arccos(np.random.uniform(size=num_pl))
	R_star = stars_chosen[:,1]
	aR = semimajor_axis_radius_ratio(P, M_star, R_star)
	b = impact_parameter(aR, alli)
    
	# Remove non-transiting planets
	itrans = (b <= 1.0)
	transP = P[itrans]
	transRp = Rp[itrans]
	transR = R_star[itrans]
	
    # Calculate SNR and detection probabilities
	dur = transit_duration(transP, aR[itrans], b[itrans])
	Ntr = num_transits(transP, Tobs[itrans], fduty[itrans])
	SNR = transit_SNR(transP, transRp, dur, transR, u1[itrans], u2[itrans], Tobs[itrans], fduty[itrans], std[itrans])
	Pdet = P_det(Ntr, SNR)
	Pwin = P_win(transP, Tobs[itrans], fduty[itrans])
    
	# Remove non-detected planets
	idet = (bernoulli.rvs(Pdet*Pwin) == 1)
	detRp = transRp[idet]
	detP = transP[idet]
	detR = transR[idet]

	# Generated observed properties
	obsdepth = transit_noise_model(detP,detRp,detR,aR[itrans][idet],b[itrans][idet],Tobs[itrans][idet],fduty[itrans][idet],std[itrans][idet])
	obsdepth[obsdepth < 0.0] = 0.0
	flag = np.random.randint(0,2,size=len(detRp))
	obsR = detR + assumed_stellar_radius(sigdown[itrans][idet],sigup[itrans][idet],flag)
	obsRp = np.sqrt(obsdepth)*obsR/rearth_in_rsol
	data_sim = np.column_stack((detP, obsRp))
	return data_sim

def distance(data_sim, params):

	'''
	This function calculates the distance between the simulated and observed data.

	Input:
		data_sim: np.array, array of simulated data
		params: dict, dictionary of parameters

	Output:
		dist: float, distance between simulated and observed data
	'''

	# Cut distances that are not in vicinity of the radius valley
	data_obs = params['dataset1']

	p_obs = data_obs[:,0]
	r_obs = data_obs[:,1]

	in_valley = (rmin <= r_obs) & (r_obs <= rmax) & (pmin <= p_obs) & (p_obs <= pmax)

	# Get period and radius 
	data_obs_valley = data_obs[in_valley]
	p_obs_valley = data_obs_valley[:,0]
	r_obs_valley = data_obs_valley[:,1]

	p_sim = data_sim[:,0]
	r_sim = data_sim[:,1]
	
	#Set up sum of occurence rate in each bin
	num_bin = num_rbins_per_run
	sum_obs = np.zeros(num_bin)
	sum_sim = np.zeros(num_bin)
	for i in range(num_bin):
		r_rearth_min, r_rearth_max = np.exp(rbins[r_start+i]), np.exp(rbins[r_start+(i+1)])
		p_days_min, p_days_max = np.exp(pbin_lower), np.exp(pbin_upper)
		sum_obs[i] = len(r_obs_valley[(p_obs_valley >= p_days_min) & (p_obs_valley <= p_days_max) & 
								(r_obs_valley >= r_rearth_min) & (r_obs_valley <= r_rearth_max)])
		sum_sim[i] = len(r_sim[(p_sim >= p_days_min) & (p_sim <= p_days_max) & 
						 (r_sim >= r_rearth_min) & (r_sim <= r_rearth_max)])
		
	#Calculate ratio of number of planets and number of stars
	sum_obs = sum_obs/Nstars
	sum_sim = sum_sim/Nstars

	#Calculate distance
	dist = []
	for i in range(num_bin):
		denominator = sum_obs[i] + sum_sim[i]
		if denominator <=1e-5:
			dist.append(1e-10)
		else:
			current_dist = np.abs(sum_obs[i]-sum_sim[i])/np.sqrt(sum_obs[i]+sum_sim[i])
			dist.append(current_dist)

	return np.atleast_1d(np.sum(dist))


    
    

