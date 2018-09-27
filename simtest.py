import multiprocessing
import os
import subprocess
import datetime

import numpy as np
from scipy.stats import bernoulli
from fileinput import filename

executable = '/home/main/chaste_build/projects/ISP/apps/HeSimulator'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

#####################
# SPSA COEFFICIENTS
#####################

a = .001
c = .1
A = 5 #10% of expected # iterations
alpha = .602
gamma = .101

###########################
# SPSA & AIC UTILITY VARS
###########################

max_iterations = 50
number_comparisons = 3000 #total number of comparison points between model output and empirical data (30 per timepoint)
number_comparisons_per_induction = number_comparisons / 3 #numberPoints must be evenly divisible by 3!
theta_spsa_size = 6; #6 parameters to optimise

#delta = vec(theta_spsaize); //vector for bernoulli random variables to perturb theta

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. lineages per loss function run
#unique sequence of RNG results for each lineage
start_seed = 0
end_seed = 99
total_seeds = end_seed - start_seed + 1
directory_name = "SPSA"
file_name = "HeSPSACounts"
log_name = "HeSPSAOutput"
deterministic_mode = 0 #1=enabled
output_mode = 0 #0=lineage counts;1=mitotic event logging;2=sequence sampling
fixture = 0 #0=He 2012;1=Wan 2016
debug_output = 0 #0=off;1=on

##########################
#GLOBAL MODEL PARAMETERS
##########################

#Values defining different marker induction timepoints & relative start time of TiL counter
earliest_lineage_start_time = 23.0 #RPCs begin to enter "He model regime" nasally at 23hpf
latest_lineage_start_time = 39.0 #last temporal retinal RPC has entered "He model regime" at 39hpf
induction_times = [ 24, 32, 48 ] 
end_time = 72.0

########################################
# SPECIFIC MODEL PARAMETERS - THETAHAT
########################################

#STOCHASTIC MITOTIC MODE
#mitotic mode per-phase probabilities
#3-phase mitotic mode time periodisation
mitotic_mode_phase_2 = 8 #These are phase lengths, so
mitotic_mode_phase_3 = 7 #Phase3 boundary = mmp2 + mmp3
phase_1_pPP = 1.0
phase_1_pPD = 0.0
phase_2_pPP = 0.2
phase_2_pPD = 0.4
phase_3_pPP = 0.2
phase_3_pPD = 0.0
he_model_params = 15

#DETERMINISTIC MITOTIC MODE
#Phase boundary shift parameters
phase_1_shape = 2
phase_1_scale = 4
phase_2_shape = 3
phase_2_scale = 4
phase_sister_shift_widths = .8
det_model_params = 12

#array "theta_spsa" is manipulated during SPSA calculations
theta_spsa = np.array([mitotic_mode_phase_2, mitotic_mode_phase_3, phase_1_pPP, phase_2_pPP, phase_2_pPD, phase_3_pPP ])

#scales ak gain sequence for probability variables
prob_scale_vector = np.array([ 1, 1, .025, .025, .025, .025 ])
shift_scale_vector = np.array([1, .5, .5, .5])

##############################
# HE ET AL EMPIRICAL RESULTS
##############################
#no. of lineages observed per induction timepoint
lineages_sampled_24 = 64
lineages_sampled_32 = 169
lineages_sampled_48 = 163

#binned histograms of counts at each timepoint
he_24_counts = np.array([ 0, 0, 1, 0, 1, 2, 0, 7, 4, 9, 6, 5, 5, 6, 5, 5, 1, 2, 2, 2, 0, 1 ])
he_32_counts = np.array([ 6, 20, 25, 22, 24, 17, 11, 7, 8, 5, 6, 3, 1, 6, 4, 1, 0, 0, 0, 2, 0, 0, 0, 1 ])
he_48_counts = np.array([ 59, 86, 2, 12, 1, 2, 0, 1 ])

#extend histogram to 1000
he_24_AIC_array = np.concatenate([he_24_counts, np.zeros(int(number_comparisons_per_induction) - he_24_counts.size)])
he_32_AIC_array = np.concatenate([he_32_counts, np.zeros(int(number_comparisons_per_induction) - he_32_counts.size)])
he_48_AIC_array = np.concatenate([he_48_counts, np.zeros(int(number_comparisons_per_induction) - he_48_counts.size)])

#Arrays containing He et al. empirical probabilities for counts 1-1000
prob_empirical_24 = np.array(he_24_AIC_array / lineages_sampled_24)
prob_empirical_32 = np.array(he_32_AIC_array / lineages_sampled_32)
prob_empirical_48 = np.array(he_48_AIC_array / lineages_sampled_48)

empirical_prob_list = [prob_empirical_24, prob_empirical_32, prob_empirical_48]

rss_plus = np.zeros(3)
rss_minus = np.zeros(3)
rss_test = np.zeros(3)

counts_plus_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts24Plus", skiprows=1, usecols=3)
counts_minus_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts24Minus", skiprows=1, usecols=3)



counts_plus_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts32Plus", skiprows=1, usecols=3)
counts_minus_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts32Minus", skiprows=1, usecols=3)


counts_plus_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts48Plus", skiprows=1, usecols=3)

counts_minus_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/SPSA/HeSPSACounts48Minus", skiprows=1, usecols=3)

counts_test_24 = np.loadtxt("/home/main/chaste_build/projects/ISP/apps/testoutput/SPSA/He24Test", skiprows=1, usecols=3)
counts_test_32 = np.loadtxt("/home/main/chaste_build/projects/ISP/apps/testoutput/SPSA/He32Test", skiprows=1, usecols=3)
counts_test_48 = np.loadtxt("/home/main/chaste_build/projects/ISP/apps/testoutput/SPSA/He48Test", skiprows=1, usecols=3)


prob_histo_plus_24, bin_edges = np.histogram(counts_plus_24, bins=1000, range=(1,1000),density=True)
prob_histo_minus_24, bin_edges = np.histogram(counts_minus_24, bins=1000, range=(1,1000),density=True)
prob_histo_plus_32, bin_edges = np.histogram(counts_plus_32, bins=1000, range=(1,1000),density=True)
prob_histo_minus_32, bin_edges = np.histogram(counts_minus_32, bins=1000, range=(1,1000),density=True)
prob_histo_plus_48, bin_edges = np.histogram(counts_plus_48, bins=1000, range=(1,1000),density=True)
prob_histo_minus_48, bin_edges = np.histogram(counts_minus_48, bins=1000, range=(1,1000),density=True)

prob_histo_test_24, bin_edges = np.histogram(counts_test_24, bins=1000, range=(1,1000),density=True)
prob_histo_test_32, bin_edges = np.histogram(counts_test_32, bins=1000, range=(1,1000),density=True)
prob_histo_test_48, bin_edges = np.histogram(counts_test_48, bins=1000, range=(1,1000),density=True)




residual_plus_24 = prob_histo_plus_24 - empirical_prob_list[0]
residual_minus_24 = prob_histo_minus_24 - empirical_prob_list[0]
residual_plus_32 = prob_histo_plus_32 - empirical_prob_list[1]
residual_minus_32 = prob_histo_minus_32 - empirical_prob_list[1]
residual_plus_48 = prob_histo_plus_48 - empirical_prob_list[2]
residual_minus_48 = prob_histo_minus_48 - empirical_prob_list[2]

residual_test_24 = prob_histo_test_24 = empirical_prob_list[0]
residual_test_32 = prob_histo_test_32 = empirical_prob_list[1]
residual_test_48 = prob_histo_test_48 = empirical_prob_list[2]

#for x in range(0,30):
#    print(residual_plus_48[x])

#for x in range(0,30):
#    print(residual_plus_24[x]-residual_test_24[x])

rss_plus[0] = np.sum(np.square(residual_plus_24))
rss_plus[1] = np.sum(np.square(residual_plus_32))
rss_plus[2] = np.sum(np.square(residual_plus_48))
rss_minus[0] = np.sum(np.square(residual_minus_24))
rss_minus[1] = np.sum(np.square(residual_minus_32))
rss_minus[2] = np.sum(np.square(residual_minus_48))

rss_test[0] = np.sum(np.square(residual_test_24))
rss_test[1] = np.sum(np.square(residual_test_32))
rss_test[2] = np.sum(np.square(residual_test_48))

#for x in range(0,len(rss_plus)):
#    print(rss_plus[x])
#    print(rss_test[x])


#print(np.sum(rss_plus))
#print(np.sum(rss_minus))
#print(np.sum(rss_test))
        
AIC_plus = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_plus))
AIC_minus = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_minus))
AIC_test = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_test))
print("AICplus " + str(AIC_plus) + " AICminus " + str(AIC_minus) + " AICtest " + str(AIC_test))
