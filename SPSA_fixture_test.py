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

a = .05
c = .5
A = 7 #10% of expected # iterations
alpha = .602
gamma = .101

###########################
# SPSA & AIC UTILITY VARS
###########################

max_iterations = 0
number_comparisons = 3000 #total number of comparison points between model output and empirical data (30 per timepoint)
number_comparisons_per_induction = number_comparisons / 3 #numberPoints must be evenly divisible by 3!

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. lineages per loss function run
#unique sequence of RNG results for each lineage
start_seed = 0
end_seed = 499
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
prob_scale_vector = np.array([ 1, 1, .05, .05, .05, .05 ])
shift_scale_vector = np.array([1, 1, 1, 1, .2])

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

#setup the log file, appending to any existing results
log_filename = "/home/main/git/chaste/projects/ISP/testoutput/" +directory_name +"/" + log_name
os.makedirs(os.path.dirname(log_filename), exist_ok=True)

log = open(log_filename,"a")
now = datetime.datetime.now()
log.write("Began SPSA optimisation of He model @\n" + str(datetime.datetime.now()) + "\n" + "k\tphase2\tphase3\tPP1\tPP2\tPD2\tPP3\n")

def main():
    global theta_spsa
    
    #SPSA algorithm iterator k starts at 0
    k=0
       
    while k <= max_iterations:
        
        #@35th iterate, increase # of seeds to decrease RNG noise
        if k == 50:
            end_seed = 4999
            total_seeds = end_seed + 1
            
        #write the parameter set to be evaluated to file
        log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4]) + "\t" + str(theta_spsa[5]) + "\n")
       
        #populate the deltak perturbation vector with samples from a .5p bernoulli +1 -1 distribution
        delta_k = bernoulli.rvs(.5, size=theta_spsa.size)
        delta_k[delta_k == 0] = -1
        
        ak = a / ((A + k + 1)**alpha) #calculate ak from gain sequence
        
        scaled_ak = ak * prob_scale_vector #scale ak appropriately for parameters expressed in hrs & percent
        ck = c / ((k + 1)**gamma) #calculate ck from gain sequence
        scaled_ck = ck * prob_scale_vector
        
        #Project theta_spsa into space bounded by ck to allow gradient sampling at bounds of probability space
        projected_theta = project_theta(theta_spsa, scaled_ck)

        #Calculate theta+ and theta- vectors for gradient estimate
        theta_plus = projected_theta + scaled_ck * delta_k
        theta_minus = projected_theta - scaled_ck * delta_k

        AIC_gradient = evaluate_AIC_gradient(k, theta_spsa, projected_theta, theta_plus, theta_minus, deterministic_mode, he_model_params, number_comparisons)

        ghat = ((AIC_gradient) / (2 * ck)) * delta_k

        log.write("ghat0: " + str(ghat[0]) + "\n")

        #update new theta_spsa
        theta_spsa = theta_spsa - scaled_ak * ghat

        #constrain updated theta_spsa
        theta_spsa = project_theta(theta_spsa, np.zeros(theta_spsa.size))

        k+=1
    
    #write final result and close log
    log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4]) + "\t" + str(theta_spsa[5]) + "\n")
    log.close

def project_theta(theta, boundary):
# Required projection is onto unit triangle modified by boundary (ck value)
# Values outside the lower bounds are reset first
# Then, for phase 2 probabilities, we project (PPa,PDa)->(PPb,PDb) such that if( (PPa+ck) + (PDa+ck) >1), PPb+ck + PDb +ck = 1.
# This takes care of the upper bound

    #project all negative parameters back into bounded space
    for i in range(0, theta.size - 1):
        if theta[i] < 0 + boundary[i]: theta[i] = boundary[i]
        i+=1

    #if PP1 or PP3 exceed 1-boundary, project back to 1-boundary
    if theta[2] > 1 - boundary[2]: theta[2] = 1 - boundary[2]
    if theta[5] > 1 - boundary[2]: theta[5] = 1 - boundary[2]
    
    
    if theta[3] + theta[4] > 1 - boundary[3]: #if pPP2 + pPD2 gives a total beyond the current boundary-reduced total of 1.0
        
        if theta[3] > (1 - 2 * boundary[2]) and theta[4] <= theta[3] - (1 - 2 * boundary[2]): # these (PP,PD) points will project into negative PD space
            theta[3] = 1 - 2 * boundary[2]
            theta[4] = boundary[2]
            
        elif theta[4] > (1 - 2 * boundary[2]) and theta[3] <= theta[4] - (1 - 2 * boundary[2]): # these (PP,PD) points will project into negative PP space
            theta[4] = 1 - 2 * boundary[2]
            theta[3] = boundary[2]
        
        else: # the (PP,PD) points can be projected onto the line PD = -PP + 1 and modified by the boundary;
            v1 = [ 1, -1 ] # vector from (0,1) to (1,0)
            v2 = [ theta[3], theta[4] - 1 ] #vector from (0,1) to (PP,PD)
            dot_product = np.dot(v1,v2)
            lengthv1 = 2
            theta[3] = (dot_product / lengthv1) - boundary[2]
            theta[4] = (1 - dot_product / lengthv1) - boundary[2]
            
    return theta

def evaluate_AIC_gradient(k, theta, projected_theta, theta_plus, theta_minus, deterministic_mode, number_parameters, total_comparisons):
    command_list = []
    base_command = executable
    
    for i in range(0,len(induction_times)):
        if deterministic_mode == 0:
            commandPlus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Plus "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(theta_plus[0])+" "\
                        +str(theta_plus[1])+" "\
                        +str(theta_plus[2])+" "\
                        +str(phase_1_pPD)+" "\
                        +str(theta_plus[3])+" "\
                        +str(theta_plus[4])+" "\
                        +str(theta_plus[5])+" "\
                        +str(phase_3_pPD)
                        
            commandMinus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Minus "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(theta_minus[0])+" "\
                        +str(theta_minus[1])+" "\
                        +str(theta_minus[2])+" "\
                        +str(phase_1_pPD)+" "\
                        +str(theta_minus[3])+" "\
                        +str(theta_minus[4])+" "\
                        +str(theta_minus[5])+" "\
                        +str(phase_3_pPD)
                        
            thetaCommand = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Theta "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(theta[0])+" "\
                        +str(theta[1])+" "\
                        +str(theta[2])+" "\
                        +str(phase_1_pPD)+" "\
                        +str(theta[3])+" "\
                        +str(theta[4])+" "\
                        +str(theta[5])+" "\
                        +str(phase_3_pPD)
            
            thetaProjected = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"ProjTheta "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(projected_theta[0])+" "\
                        +str(projected_theta[1])+" "\
                        +str(projected_theta[2])+" "\
                        +str(phase_1_pPD)+" "\
                        +str(projected_theta[3])+" "\
                        +str(projected_theta[4])+" "\
                        +str(projected_theta[5])+" "\
                        +str(phase_3_pPD)
                        
            command_list.append(commandPlus)
            command_list.append(commandMinus)
            command_list.append(thetaCommand)
            command_list.append(thetaProjected)

        if deterministic_mode == 1:
            commandPlus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Plus "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(theta_plus[0])+" "\
                        +str(theta_plus[1])+" "\
                        +str(theta_plus[2])+" "\
                        +str(theta_plus[3])+" "\
                        +str(theta_plus[4])
                        
            commandMinus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Minus "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "\
                        +str(theta_minus[0])+" "\
                        +str(theta_minus[1])+" "\
                        +str(theta_minus[2])+" "\
                        +str(theta_minus[3])+" "\
                        +str(theta_minus[4])
                        
            command_list.append(commandPlus)
            command_list.append(commandMinus)
            
    
    # Use processes equal to the number of cpus available
    cpu_count = multiprocessing.cpu_count()

    print("Starting simulations for iterate " + str(k) + " with " + str(cpu_count) + " processes")

    log.flush() #required to prevent pool from jamming up log for some reason
    
    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=cpu_count)

    # Pass the list of bash commands to the pool, block until pool is complete
    pool.map(execute_command, command_list, 1)
    
    rss_plus = np.zeros(len(induction_times))
    rss_minus = np.zeros(len(induction_times))
    rss_theta = np.zeros(len(induction_times))
    rss_thetap = np.zeros(len(induction_times))

    
    for i in range(0,len(induction_times)):
        counts_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Plus", skiprows=1, usecols=3)
        counts_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Minus", skiprows=1, usecols=3)
        counts_theta = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Theta", skiprows=1, usecols=3)
        counts_thetap =  np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "ProjTheta", skiprows=1, usecols=3)
        prob_histo_plus, bin_edges = np.histogram(counts_plus, bins=1000, range=(1,1000),density=True)
        prob_histo_minus, bin_edges = np.histogram(counts_minus, bins=1000, range=(1,1000),density=True)
        prob_theta, bin_edges = np.histogram(counts_theta, bins=1000, range=(1,1000),density=True)
        prob_thetap, bin_edges = np.histogram(counts_thetap, bins=1000, range=(1,1000),density=True)
        residual_plus = prob_histo_plus - empirical_prob_list[i]
        residual_minus = prob_histo_minus - empirical_prob_list[i]
        residual_theta = prob_theta - empirical_prob_list[i]
        residual_thetap = prob_thetap - empirical_prob_list[i]
        rss_plus[i] = np.sum(np.square(residual_plus))
        rss_minus[i] = np.sum(np.square(residual_minus))
        rss_theta[i] = np.sum(np.square(residual_theta))
        rss_thetap[i] = np.sum(np.square(residual_thetap))
        
    AIC_plus = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_plus))
    AIC_minus = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_minus))
    AIC_theta = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_theta))
    AIC_thetap = 2 * he_model_params + number_comparisons * np.log(np.sum(rss_thetap))

    log.write("theta_plus:\n")
    log.write(str(k) + "\t" + str(theta_plus[0]) + "\t" + str(theta_plus[1]) + "\t" + str(theta_plus[2]) + "\t" + str(theta_plus[3]) + "\t" + str(theta_plus[4]) + "\t" + str(theta_plus[5]) + "\n")
    log.write("theta_minus:\n")
    log.write(str(k) + "\t" + str(theta_minus[0]) + "\t" + str(theta_minus[1]) + "\t" + str(theta_minus[2]) + "\t" + str(theta_minus[3]) + "\t" + str(theta_minus[4]) + "\t" + str(theta_minus[5]) + "\n")
    log.write("ThetaAIC: " + str(AIC_theta) + " ProjThetaAIC: "+str(AIC_thetap)+" PositiveAIC: " + str(AIC_plus) + " NegativeAIC: " + str(AIC_minus) + "\n")

    AIC_gradient_sample = AIC_plus - AIC_minus;

    return AIC_gradient_sample

# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
