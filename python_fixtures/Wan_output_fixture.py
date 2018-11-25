import multiprocessing
import os
import subprocess
import datetime

import numpy as np
from openpyxl.styles.builtins import output
from sklearn.datasets.tests.test_svmlight_format import currdir

executable = '/home/main/chaste_build/projects/ISP/apps/WanSimulator'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. simulated CMZs per run
#unique sequence of RNG results for each lineage
start_seed = 0
end_seed = 99 #total seeds should be divisible by # cores

output_directory = "WanOutput"

##########################
#MODEL PARAMETERS
##########################

#CMZ residency time 
wan_residency_time = 17.0
#factor to divide 3dpf progenitor population by to obtain estimate of stem cell population
stem_divisor = 10
#3dpf CMZ progenitor population mean and standard deviation
progenitor_mean = 792
progenitor_std = 160

#stem cell cycle parameters- cell cycle RV is a shifted gamma distribution
stem_gamma_shift = 4
stem_gamma_shape = 6.5 #mean 30 hr stem cell time - in reality is probably more like 60+
stem_gamma_scale = 4

progenitor_gamma_shift = 4  #default He et al. values, mean 6 hr cycle time
progenitor_gamma_shape = 2
progenitor_gamma_scale = 1
progenitor_gamma_sister = 1

cmz_theta_string = str(wan_residency_time) + " " + str(stem_divisor) + " " + str(progenitor_mean) + " " + str(progenitor_std) + " " + str(stem_gamma_shift) + " " + str(stem_gamma_shape) + " " + str(stem_gamma_scale) + " " + str(progenitor_gamma_shift) + " " + str(progenitor_gamma_shape) + " " + str(progenitor_gamma_scale) + " " + str(progenitor_gamma_sister)

#STOCHASTIC MITOTIC MODE
#HE ORIGINAL PARAMETERS
#mitotic mode per-phase probabilities
#3-phase mitotic mode time periodisation
mitotic_mode_phase_2 = 8 #These are phase boundaries rather than lengths as in eg. He_output_fixture.py
mitotic_mode_phase_3 = 15
phase_1_pPP = 1.0
phase_1_pPD = 0.0
phase_2_pPP = 0.2
phase_2_pPD = 0.4
phase_3_pPP = 0.2
phase_3_pPD = 0.0

stochastic_theta_string = str(mitotic_mode_phase_2) + " " + str(mitotic_mode_phase_3) + " " +str(phase_1_pPP) + " " + str(phase_1_pPD) + " " + str(phase_2_pPP) + " " + str(phase_2_pPD) + " " + str(phase_3_pPP) + " " + str(phase_3_pPD)

def main():

    command_list = []

    # Use processes equal to the number of cpus available
    cpu_count = multiprocessing.cpu_count()
    seeds_per_command = (end_seed + 1) / cpu_count
    base_command = executable + " " + output_directory
    
    for i in range(0,cpu_count):
        curr_start_seed = i * seeds_per_command
        curr_end_seed = curr_start_seed + (seeds_per_command - 1)
        
        command = base_command + " "\
            +str(curr_start_seed) + " "\
            +str(curr_end_seed) + " "\
            +cmz_theta_string + " "\
            +stochastic_theta_string
        
        #command_list.append(command)
        
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 18.0 24.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 42.0 49.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 51.0 57.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 94.0 99.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 58.0 64.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 65.0 70.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")
    command_list.append("/home/main/chaste_build/projects/ISP/apps/WanSimulator WanOutput 71.0 74.0 17.0 10 792 160 4 6.5 4 4 2 1 1 8 15 1.0 0.0 0.2 0.4 0.2 0.0")

    print("Starting simulations with " + str(cpu_count) + " processes")
   
    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=cpu_count)

    # Pass the list of bash commands to the pool, block until pool is complete
    pool.map(execute_command, command_list, 1)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    print("Executing command: " + cmd)
    return subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
