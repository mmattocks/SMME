import multiprocessing
import os
import subprocess
import datetime

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli
from fileinput import filename
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

executable = '/home/main/chaste_build/projects/ISP/apps/HeSimulator'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

#####################
# SPSA COEFFICIENTS
#####################

a_s = .025
c_s = .1

a_d = .0019
c_d = .1

A = 20 #10% of expected # iterations
alpha = .602
gamma = .101

###########################
# SPSA & AIC UTILITY VARS
###########################

max_iterations = 199
number_comparisons = 3030 #total number of comparison points between model output and empirical data (1000 per induction time, 10 per rate mode)
number_comparisons_per_induction = 1000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. lineages per loss function run
#unique sequence of RNG results for each lineage
start_seed = 0
end_seed = 249
rate_end_seed = 99
directory_name = "SPSA"
file_name_he = "HeSPSA"
file_name_det = "DetSPSA"
log_name = "HeSPSAOutput"
deterministic_modes = [0, 1] #1=det. mode enabled
count_output_mode = 0 #0=lineage counts;1=mitotic event logging;2=sequence sampling
event_output_mode = 1
fixture = 0 #0=He 2012;1=Wan 2016
debug_output = 0 #0=off;1=on
ath5founder = 0 #0=no morpholino 1=ath5 morpholino

##########################
#GLOBAL MODEL PARAMETERS
##########################

#Values defining different marker induction timepoints & relative start time of TiL counter
earliest_lineage_start_time = 23.0 #RPCs begin to enter "He model regime" nasally at 23hpf
latest_lineage_start_time = 39.0 #last temporal retinal RPC has entered "He model regime" at 39hpf
induction_times = [ 24, 32, 48 ] 
end_time = 72.0
rate_end_time = 80

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
phase_1_shape = 3
phase_1_scale = 2
phase_2_shape = 2
phase_2_scale = 2
phase_sister_shift_widths = .25
phase_offset = 0
det_model_params = 13

#array "theta_spsa" is manipulated during SPSA calculations, these are starting point references
stochastic_theta_zero = np.array([mitotic_mode_phase_2, mitotic_mode_phase_3, phase_2_pPP, phase_2_pPD, phase_3_pPP ])
determinsitic_theta_zero = np.array([phase_1_shape, phase_1_scale, phase_2_shape, phase_2_scale, phase_sister_shift_widths, phase_offset])

#scales ak gain sequence for probability variables
prob_scale_vector = np.array([ 1, 1, .1, .1, .1])
shift_scale_vector = np.array([1, 1, 1, 1, 1, 3])

##############################
# HE ET AL EMPIRICAL RESULTS
##############################

count_bin_sequence = np.arange(1,32,1)
count_x_sequence = np.arange(1,31,1)
count_trim_value = 30

rate_bin_sequence = np.arange(30,85,5)
rate_x_sequence = np.arange(30,80,5)
rate_trim_value = 10

#Count probability arrays
raw_counts = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_counts.csv', skiprows=1, usecols=(3,4,5,6,7,8,9,10)) #collect the per-cell-type counts
raw_counts_24 = raw_counts[0:64,:]
raw_counts_32 = raw_counts[64:233,:]
raw_counts_48 = raw_counts[233:396,:]

prob_24,bin_edges = np.histogram(np.sum(raw_counts_24,axis=1),count_bin_sequence,density=True) #obtain probability density histogram for counts, retabulating by summing across types
prob_32,bin_edges = np.histogram(np.sum(raw_counts_32,axis=1),count_bin_sequence,density=True)
prob_48,bin_edges = np.histogram(np.sum(raw_counts_48,axis=1),count_bin_sequence,density=True)

#extend histogram to 1000 - allows large lineages of some iterates to be included in AIC comparison
prob_empirical_24 = np.concatenate([prob_24, np.zeros(int(number_comparisons_per_induction) - prob_24.size)])
prob_empirical_32 = np.concatenate([prob_32, np.zeros(int(number_comparisons_per_induction) - prob_32.size)])
prob_empirical_48 = np.concatenate([prob_48, np.zeros(int(number_comparisons_per_induction) - prob_48.size)])

count_prob_list = [prob_empirical_24, prob_empirical_32, prob_empirical_48]

#no. of lineages observed per induction timepoint/event group
lineages_sampled_24 = 64
lineages_sampled_32 = 169
lineages_sampled_48 = 163
lineages_sampled_events = 60
lineages_sampled_list = [lineages_sampled_24,lineages_sampled_32,lineages_sampled_48]

#Mitotic mode rate probability arrays
raw_events = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_lineages.csv', skiprows=1, usecols=(3,5,8))
observed_events = raw_events[np.where(raw_events[:,2]==1)] #exclude any mitosis whose time was too early for recording

observed_PP = observed_events[np.where(observed_events[:,0]==0)]
observed_PD = observed_events[np.where(observed_events[:,0]==1)]
observed_DD = observed_events[np.where(observed_events[:,0]==2)]

histo_PP,bin_edges = np.histogram(observed_PP,rate_bin_sequence,density=False)
histo_PD,bin_edges = np.histogram(observed_PD,rate_bin_sequence,density=False)
histo_DD,bin_edges = np.histogram(observed_DD,rate_bin_sequence,density=False)

#hourly per-lineage probabilities- NOT probability density function
prob_empirical_PP = np.array((histo_PP/lineages_sampled_events)/5)
prob_empirical_PD = np.array((histo_PD/lineages_sampled_events)/5)
prob_empirical_DD = np.array((histo_DD/lineages_sampled_events)/5)

event_prob_list = [prob_empirical_PP,prob_empirical_PD,prob_empirical_DD]

#setup the log file, appending to any existing results
log_filename = "/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" +directory_name +"/" + log_name
os.makedirs(os.path.dirname(log_filename), exist_ok=True)

log = open(log_filename,"w")

#plotting utility stuff
#interactive plot mode
plt.ion()
#setup new iterate plot
fig, ((plt0, plt1, plt2),(plt3, plt4, plt5)) = plt.subplots(2,3,figsize=(12,6))
plot_list = [plt0,plt1,plt2,plt3,plt4,plt5]


def main():
    global theta_spsa, a,c,file_name,scale_vector,end_seed,rate_end_seed,alpha,gamma
    
    for m in range(0,len(deterministic_modes)):
        
        #traceable RNG
        p_RNG = np.random.RandomState(seed=786)

        deterministic_mode = deterministic_modes[m]
        if deterministic_mode == 0:
            theta_spsa = stochastic_theta_zero
            a = a_s
            c = c_s
            file_name = file_name_he
            scale_vector = prob_scale_vector
            now = datetime.datetime.now()
            log.write("Began SPSA optimisation of He model @\n" + str(datetime.datetime.now()) + "\n")
            log.write("k\tphase2\tphase3\tPP2\tPD2\tPP3\n")
            number_params = he_model_params
        if deterministic_mode == 1:
            theta_spsa = determinsitic_theta_zero
            a = a_d
            c = c_d
            file_name = file_name_det
            scale_vector = shift_scale_vector
            now = datetime.datetime.now()
            log.write("Began SPSA optimisation of deterministic model @\n" + str(datetime.datetime.now()) + "\n")
            log.write("k\tp1Sh\tp1Sc\tp2Sh\tp2Sc\tsisterShift\toffset\n")
            number_params = det_model_params
        
        #SPSA algorithm iterator k starts at 0
        k=0
       
        while k <= max_iterations:
            if k == 0:
                end_seed = 249
                rate_end_seed = 99
        
            #@defined iterate, increase # of seed to decrease RNG noise to low level
            if k == 169:
                end_seed = 999
                rate_end_seed = 249
        
        
            #@defined iterate, increase # of seeds to decrease RNG noise to close to nil, switch to asymptotically optimal alpha and gamma
            if k == 189:
                end_seed = 4999
                rate_end_seed = 1249
                alpha = 1
                gamma = (1/6)
            
            #write the parameter set to be evaluated to file
            if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4]) + "\n")
            if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4])+ "\t" + str(theta_spsa[5]) + "\n")
        
            #populate the deltak perturbation vector with samples from a .5p bernoulli +1 -1 distribution
            delta_k = bernoulli.rvs(.5, size=theta_spsa.size, random_state=p_RNG)
            delta_k[delta_k == 0] = -1
        
            ak = a / ((A + k + 1)**alpha) #calculate ak from gain sequence
            scaled_ak = ak * scale_vector #scale ak appropriately for parameters expressed in hrs & percent

            ck = c / ((k + 1)**gamma) #calculate ck from gain sequence
            scaled_ck = ck * scale_vector
        
            #Project theta_spsa into space bounded by ck to allow gradient sampling at bounds of probability space
            projected_theta = project_theta(theta_spsa, scaled_ck, deterministic_mode)

            #Calculate theta+ and theta- vectors for gradient estimate
            theta_plus = projected_theta + scaled_ck * delta_k
            theta_minus = projected_theta - scaled_ck * delta_k

            AIC_gradient = evaluate_AIC_gradient(k, theta_plus, theta_minus, deterministic_mode, number_params, file_name)

            ghat = ((AIC_gradient) / (2 * ck)) * delta_k

            log.write("ghat0: " + str(ghat[0]) + "\n")

            #update new theta_spsa
            theta_spsa = theta_spsa - scaled_ak * ghat

            #constrain updated theta_spsa
            theta_spsa = project_theta(theta_spsa, np.zeros(theta_spsa.size), deterministic_mode)

            k+=1
        
        #write final result and close log
        if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4]) + "\n")
        if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_spsa[0]) + "\t" + str(theta_spsa[1]) + "\t" + str(theta_spsa[2]) + "\t" + str(theta_spsa[3]) + "\t" + str(theta_spsa[4])+ "\t" + str(theta_spsa[5]) + "\n")

    log.close

def project_theta(theta, boundary, deterministic_mode):
# Required projection is onto unit triangle modified by boundary (ck value)
# Values outside the lower bounds are reset first
# Then, for phase 2 probabilities, we project (PPa,PDa)->(PPb,PDb) such that if( (PPa+ck) + (PDa+ck) >1), PPb+ck + PDb +ck = 1.
# This takes care of the upper bound

    #project all negative parameters back into bounded space
    if deterministic_mode == False:
        for i in range(0, theta.size):
            if i == 0:
                if theta[i] < boundary[i] + 4.0: theta[i] = boundary[i] + 4.0 #project mitotic_mode_phase_2 to a minimum of 4- below this param has no effect (due to refractory period after first division)
            else:
                if theta[i] < boundary[i]: theta[i] = boundary[i]
                
        #if PP3 exceeds 1-boundary, project back to 1-boundary
        if theta[4] > 1 - boundary[4]: theta[4] = 1 - boundary[4]
        
        if theta[2] + theta[3] > 1 - boundary[2]: #if pPP2 + pPD2 gives a total beyond the current boundary-reduced total of 1.0
        
            if theta[2] > (1 - 2 * boundary[2]) and theta[2] <= theta[3] - (1 - 2 * boundary[2]): # these (PP,PD) points will project into negative PD space
                theta[2] = 1 - 2 * boundary[2]
                theta[3] = boundary[2]
            
            elif theta[3] > (1 - 2 * boundary[2]) and theta[2] <= theta[3] - (1 - 2 * boundary[2]): # these (PP,PD) points will project into negative PP space
                theta[3] = 1 - 2 * boundary[2]
                theta[2] = boundary[2]
        
            else: # the (PP,PD) points can be projected onto the line PD = -PP + 1 and modified by the boundary;
                v1 = [ 1, -1 ] # vector from (0,1) to (1,0)
                v2 = [ theta[2], theta[3] - 1 ] #vector from (0,1) to (PP,PD)
                dot_product = np.dot(v1,v2)
                lengthv1 = 2
                theta[2] = (dot_product / lengthv1) - boundary[2]
                theta[3] = (1 - dot_product / lengthv1) - boundary[2]

    if deterministic_mode == True:
        for i in range(0, theta.size-1): #exclude offset param
            if theta[i] < boundary[i]: theta[i] = boundary[i] + 0.1 #prevents non->0 arguments for deterministic mode

        #projects sister shift value such that 95% of sister shift values will be less than smallest mean phase time
        minPhase = min(theta[0]*theta[1],theta[2]*theta[3])
        if theta[4] > minPhase/2 - boundary[4]: theta[4] = minPhase/2 - boundary[4]

    return theta

def evaluate_AIC_gradient(k, theta_plus, theta_minus, deterministic_mode, number_params, file_name):
    #Form the simulator commands for current thetas and deterministic modes
    command_list = []
    base_command = executable
    
    rate_settings = str(event_output_mode)+" "\
                    +str(deterministic_mode)+" "\
                    +str(fixture)+" "\
                    +str(ath5founder)+" "\
                    +str(debug_output)+" "\
                    +str(start_seed)+" "\
                    +str(rate_end_seed)+" "\
                    +str(earliest_lineage_start_time)+" "\
                    +str(earliest_lineage_start_time)+" "\
                    +str(latest_lineage_start_time)+" "\
                    +str(rate_end_time)+" "
                    
    count_settings_1 = str(count_output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(start_seed)+" "\
                        +str(end_seed)+" "
    
    count_settings_2 = str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(end_time)+" "
    
    if deterministic_mode == 0:
        
        stochastic_params_plus = str(theta_plus[0])+" "\
                    +str(theta_plus[1])+" "\
                    +str(phase_1_pPP)+" "\
                    +str(phase_1_pPD)+" "\
                    +str(theta_plus[2])+" "\
                    +str(theta_plus[3])+" "\
                    +str(theta_plus[4])+" "\
                    +str(phase_3_pPD)
        
        stochastic_params_minus = str(theta_minus[0])+" "\
                    +str(theta_minus[1])+" "\
                    +str(phase_1_pPP)+" "\
                    +str(phase_1_pPD)+" "\
                    +str(theta_minus[2])+" "\
                    +str(theta_minus[3])+" "\
                    +str(theta_minus[4])+" "\
                    +str(phase_3_pPD)
        
        command_rate_plus = base_command\
                    +" "+directory_name+" "+ file_name+"RatePlus "\
                    +rate_settings\
                    +stochastic_params_plus
                        
        command_rate_minus = base_command\
                    +" "+directory_name+" "+ file_name+"RateMinus "\
                    +rate_settings\
                    +stochastic_params_minus

        command_list.append(command_rate_plus)
        command_list.append(command_rate_minus)
        
        for i in range(0,len(induction_times)):
            command_plus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Plus "\
                        +count_settings_1\
                        +str(induction_times[i])+" "\
                        +count_settings_2\
                        +stochastic_params_plus
                        
            command_minus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Minus "\
                        +count_settings_1\
                        +str(induction_times[i])+" "\
                        +count_settings_2\
                        +stochastic_params_minus
                        
            command_list.append(command_plus)
            command_list.append(command_minus)
                        
    if deterministic_mode == 1:
        
        deterministic_params_plus = str(theta_plus[0])+" "\
                    +str(theta_plus[1])+" "\
                    +str(theta_plus[2])+" "\
                    +str(theta_plus[3])+" "\
                    +str(theta_plus[4])+" "\
                    +str(theta_plus[5])
                    
        deterministic_params_minus = str(theta_minus[0])+" "\
                    +str(theta_minus[1])+" "\
                    +str(theta_minus[2])+" "\
                    +str(theta_minus[3])+" "\
                    +str(theta_minus[4])+" "\
                    +str(theta_minus[5])
        
        command_rate_plus = base_command\
                    +" "+directory_name+" "+ file_name+"RatePlus "\
                    +rate_settings\
                    +deterministic_params_plus
                
        command_rate_minus = base_command\
                    +" "+directory_name+" "+ file_name+"RateMinus "\
                    +rate_settings\
                    +deterministic_params_minus

        command_list.append(command_rate_plus)
        command_list.append(command_rate_minus)
        
        for i in range(0,len(induction_times)):
            command_plus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Plus "\
                        +count_settings_1\
                        +str(induction_times[i])+" "\
                        +count_settings_2\
                        +deterministic_params_plus
                        
            command_minus = base_command\
                        +" "+directory_name+" "+file_name+str(induction_times[i])+"Minus "\
                        +count_settings_1\
                        +str(induction_times[i])+" "\
                        +count_settings_2\
                        +deterministic_params_minus
                        
            command_list.append(command_plus)
            command_list.append(command_minus)
            

    # Use processes equal to the number of cpus available
    cpu_count = multiprocessing.cpu_count()

    print("Starting simulations for iterate " + str(k) + " with " + str(cpu_count) + " processes, " + str(end_seed+1) + " lineages simulated for counts, " + str(rate_end_seed+1) + " lineages for events")

    log.flush() #required to prevent pool from jamming up log for some reason
    
    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=cpu_count)

    # Pass the list of bash commands to the pool, block until pool is complete
    pool.map(execute_command, command_list, 1)
    
    #these numpy arrays hold the individual RSS vals for timepoints + rates
    rss_plus = np.zeros(len(induction_times)+3)
    rss_minus = np.zeros(len(induction_times)+3)
    
    #clear the plots on the interactive figure
    for i in range(0, len(plot_list)):
        plt.sca(plot_list[i])
        plt.cla()
    
    #calculate RSS for each induction timepoint
    for i in range(0,len(induction_times)):
        counts_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Plus", skiprows=1, usecols=3)
        counts_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Minus", skiprows=1, usecols=3)
        prob_histo_plus, bin_edges = np.histogram(counts_plus, bins=1000, range=(1,1001),density=True)
        prob_histo_minus, bin_edges = np.histogram(counts_minus, bins=1000, range=(1,1001),density=True)
        
        residual_plus = prob_histo_plus - count_prob_list[i]
        residual_minus = prob_histo_minus - count_prob_list[i]
        
        rss_plus[i] = np.sum(np.square(residual_plus))
        rss_minus[i] = np.sum(np.square(residual_minus))
        
        plotter(plot_list[i], counts_plus, counts_minus, count_prob_list[i],lineages_sampled_list[i], 0)

        
    rates_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" + directory_name + "/" + file_name + "RatePlus", skiprows=1, usecols=(0,3))
    rates_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" + directory_name + "/" + file_name + "RateMinus", skiprows=1, usecols=(0,3))
    
    for i in range (0,3):
        mode_rate_plus = np.array(rates_plus[np.where(rates_plus[:,1]==i)])
        mode_rate_minus = np.array(rates_minus[np.where(rates_minus[:,1]==i)])

        histo_rate_plus, bin_edges = np.histogram(mode_rate_plus[:,0], rate_bin_sequence, density=False)
        histo_rate_minus, bin_edges = np.histogram(mode_rate_minus[:,0], rate_bin_sequence, density=False)
        
        #hourly per-lineage probabilities
        prob_histo_rate_plus = np.array(histo_rate_plus / ((rate_end_seed + 1)*5))
        prob_histo_rate_minus = np.array(histo_rate_plus / ((rate_end_seed + 1)*5))
        
        residual_plus = prob_histo_rate_plus - event_prob_list[i]
        residual_minus = prob_histo_rate_minus - event_prob_list[i]
        
        #PD RESIUDAL WEIGHTING
        if i == 1:
            rss_plus[i+3] = np.sum(1.5*np.square(residual_plus))
            rss_minus[i+3] = np.sum(1.5*np.square(residual_minus))
        
        else:
            rss_plus[i+3] = np.sum(np.square(residual_plus))
            rss_minus[i+3] = np.sum(np.square(residual_minus))
        
        plotter(plot_list[i+3], mode_rate_plus[:,0], mode_rate_minus[:,0], event_prob_list[i],lineages_sampled_events, 1)

    AIC_plus = 2 * number_params + number_comparisons * np.log(np.sum(rss_plus))
    AIC_minus = 2 * number_params + number_comparisons * np.log(np.sum(rss_minus))
    
    plt.show()
    plt.savefig("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/" +directory_name +"/" + "SDmode" + str(deterministic_mode) + "iterate" + str(k) + ".png")
    
    log.write("theta_plus:\n")
    if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_plus[0]) + "\t" + str(theta_plus[1]) + "\t" + str(theta_plus[2]) + "\t" + str(theta_plus[3]) + "\t" + str(theta_plus[4]) + "\n")
    if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_plus[0]) + "\t" + str(theta_plus[1]) + "\t" + str(theta_plus[2]) + "\t" + str(theta_plus[3]) + "\t" + str(theta_plus[4]) + "\t" + str(theta_plus[5]) +"\n")
    log.write("theta_minus:\n")
    if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_minus[0]) + "\t" + str(theta_minus[1]) + "\t" + str(theta_minus[2]) + "\t" + str(theta_minus[3]) + "\t" + str(theta_minus[4]) + "\n")
    if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_minus[0]) + "\t" + str(theta_minus[1]) + "\t" + str(theta_minus[2]) + "\t" + str(theta_minus[3]) + "\t" + str(theta_minus[4]) + "\t" + str(theta_minus[5]) + "\n")

    log.write("PositiveAIC: " + str(AIC_plus) + " NegativeAIC: " + str(AIC_minus) + "\n")

    AIC_gradient_sample = AIC_plus - AIC_minus;

    return AIC_gradient_sample

#data plotter function for monitoring SPSA results
def plotter(subplot,plus,minus,empirical_prob,samples,mode):
    global plt0,plt1,plt2,plt3,plt4,plt5, prob_histo_plus, prob_histo_minus
    bin_sequence = []
    x_sequence = []
    trim_value = 0
    
    if mode == 0:
        bin_sequence = count_bin_sequence
        x_sequence = count_x_sequence
        trim_value = count_trim_value
        prob_histo_plus, bin_edges = np.histogram(plus, bin_sequence, density=True)
        prob_histo_minus, bin_edges = np.histogram(minus, bin_sequence, density=True)
    
    if mode == 1:
        bin_sequence = rate_bin_sequence
        x_sequence = rate_x_sequence
        trim_value = rate_trim_value
        histo_plus, bin_edges = np.histogram(plus, bin_sequence, density=False)
        histo_minus, bin_edges = np.histogram(minus, bin_sequence, density=False)
        prob_histo_plus = np.array(histo_plus / ((rate_end_seed + 1)*5))
        prob_histo_minus = np.array(histo_plus / ((rate_end_seed + 1)*5))
    
    trimmed_prob = empirical_prob[0:trim_value]
    
    subplot.plot(x_sequence,trimmed_prob, 'k+')
    plt.pause(0.0001)
    
    interval_plus = sampler(plus,samples,bin_sequence)
    
    subplot.plot(x_sequence,prob_histo_plus, 'g-')
    plt.pause(0.0001)
    subplot.fill_between(x_sequence, (prob_histo_plus - interval_plus), (prob_histo_plus + interval_plus), alpha=0.2, edgecolor='#008000', facecolor='#00FF00')
    plt.pause(0.0001)
    
    interval_minus = sampler(minus,samples,bin_sequence)
    
    subplot.plot(x_sequence,prob_histo_minus, 'm-')
    plt.pause(0.0001)
    subplot.fill_between(x_sequence, (prob_histo_minus - interval_minus), (prob_histo_minus + interval_minus), alpha=0.2, edgecolor='#800080', facecolor='#FF00FF')
    plt.pause(0.0001)
    
    if subplot==plt0: plt0.set_ylim((0,.20))
    if subplot==plt1: plt1.set_ylim(0,.20)
    if subplot==plt2: plt2.set_ylim(0,.7)
    if subplot==plt3: plt3.set_ylim(0, .30)
    if subplot==plt4: plt4.set_ylim(0, .15)
    if subplot==plt5: plt5.set_ylim(0, .35)
    
def sampler(data,samples,bin_sequence):
    
    if len(data) == 0: #catches edge case w/ no entries for a mitotic mode
        data=[0]
    
    base_sample=np.zeros((error_samples,len(bin_sequence)-1))

    for i in range(0,error_samples):
        new_data_sample = np.random.choice(data,samples)
        new_histo_prob, bin_edges = np.histogram(new_data_sample, bin_sequence, density=True)
        base_sample[i,:] = new_histo_prob
        
    sample_95CI = np.array(2 * (np.std(base_sample,0)))
    
    return sample_95CI

# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
