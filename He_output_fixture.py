import multiprocessing
import os
import subprocess
import datetime

import numpy as np
from imageio.plugins._bsdf import BsdfSerializer

executable = '/home/main/chaste_build/projects/ISP/apps/HeSimulator'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. lineages per loss function run
#unique sequence of RNG results for each lineage
start_seed = 0
end_seed_counts = 99999
end_seed_events = 9999
deterministic_modes = [1,0] #1=enabled
debug_output = 0 #0=off;1=on

count_directory = "HeCounts"
count_filenames = ["induction", "wan", "ath5", "validate"]
event_directory = "HeModeEvent"
event_filenames = ["inductionMode", "validateMode"]
induction_times = [24, 32, 48]

##########################
#GLOBAL MODEL PARAMETERS
##########################

#Values defining different marker induction timepoints & relative start time of TiL counter
earliest_lineage_start_time = 23.0 #RPCs begin to enter "He model regime" nasally at 23hpf
latest_lineage_start_time = 39.0 #last temporal retinal RPC has entered "He model regime" at 39hpf
counts_end_time = 72.0
events_end_time = 80.0
wan_end_time = 17.0

########################################
# SPECIFIC MODEL PARAMETERS - THETAHAT
########################################

#These parameters are the results of the SPSA optimisation fixture

#STOCHASTIC MITOTIC MODE
#mitotic mode per-phase probabilities
#3-phase mitotic mode time periodisation
mitotic_mode_phase_2 = 2.997 #These are phase lengths, so
mitotic_mode_phase_3 = 9.730 #Phase3 boundary = mmp2 + mmp3
phase_1_pPP = 1.0
phase_1_pPD = 0.0
phase_2_pPP = 0.2728
phase_2_pPD = 0.4859
phase_3_pPP = 0.2843
phase_3_pPD = 0.0
he_model_params = 15

#DETERMINISTIC MITOTIC MODE
#Phase boundary shift parameters
phase_1_shape = 2.789
phase_1_scale = 1.430
phase_2_shape = 3.830
phase_2_scale = 3.524
phase_sister_shift_widths = .7288
phase_offset = .07611
det_model_params = 12

stochastic_theta_string = str(mitotic_mode_phase_2) + " " + str(mitotic_mode_phase_3) + " " +str(phase_1_pPP) + " " + str(phase_1_pPD) + " " + str(phase_2_pPP) + " " + str(phase_2_pPD) + " " + str(phase_3_pPP) + " " + str(phase_3_pPD)
deterministic_theta_string = str(phase_1_shape) + " " + str(phase_1_scale) + " " + str(phase_2_shape) + " " + str(phase_2_scale) + " " + str(phase_sister_shift_widths) + " " + str(phase_offset)

def main():

    command_list = []

    for m in range(0,len(deterministic_modes)):
        curr_list = []
        
        deterministic_mode = deterministic_modes[m]
        command_count = 0 #for iterating seed numbers
        base_command = executable
        
        #induction count commands
        for i in range(0,len(induction_times)):
            
            output_mode = 0 #0=lineage counts;1=mitotic event logging;2=sequence sampling
            fixture = 0 #0=He 2012;1=Wan 2016
            curr_start_seed = start_seed + command_count * (end_seed_counts+1)
            curr_end_seed = curr_start_seed + end_seed_counts
            ath5founder = 0
            
            command = base_command+" "\
                        +count_directory+" "\
                        +count_filenames[0]+str(induction_times[i])+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(induction_times[i])+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(counts_end_time)+" "
            curr_list.append(command)
            command_count += 1
        
        #wan command
        output_mode = 0
        fixture = 1 
        curr_start_seed = start_seed + command_count * (end_seed_counts+1)
        curr_end_seed = curr_start_seed + end_seed_counts
        ath5founder = 0
        wan_command = base_command+" "\
                        +count_directory+" "\
                        +count_filenames[1]+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(0)+" "\
                        +str(0)+" "\
                        +str(1)+" "\
                        +str(wan_end_time)+" "
        curr_list.append(wan_command)
        command_count += 1
            
        #ath5 command
        output_mode = 0
        fixture = 2 
        curr_start_seed = start_seed + command_count * (end_seed_counts+1)
        curr_end_seed = curr_start_seed + end_seed_counts
        ath5founder = 1

        ath5_command = base_command+" "\
                        +count_directory+" "\
                        +count_filenames[2]+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(0)+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(counts_end_time)+" "
        curr_list.append(ath5_command)
        command_count += 1
        
        #validate counts command
        curr_start_seed = start_seed + command_count * (end_seed_counts+1)
        curr_end_seed = curr_start_seed + end_seed_counts
        ath5founder = 0

        validate_command = base_command+" "\
                        +count_directory+" "\
                        +count_filenames[3]+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(0)+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(counts_end_time)+" "
        curr_list.append(validate_command)
        command_count += 1
        
        #mitotic mode rate command
        output_mode = 1
        fixture = 0
        curr_start_seed = start_seed + command_count * (end_seed_counts+1)
        curr_end_seed = curr_start_seed + end_seed_events
        mode_rate_command = base_command+" "\
                        +event_directory+" "\
                        +event_filenames[0]+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(events_end_time)+" "
        
        curr_list.append(mode_rate_command)
        command_count += 1

        #validate rate command
        fixture = 2
        curr_start_seed = start_seed + command_count * (end_seed_counts+1)
        curr_end_seed = curr_start_seed + end_seed_events
        validate_rate_command = base_command+" "\
                        +event_directory+" "\
                        +event_filenames[1]+"SDMode"+str(deterministic_mode)+" "\
                        +str(output_mode)+" "\
                        +str(deterministic_mode)+" "\
                        +str(fixture)+" "\
                        +str(ath5founder)+" "\
                        +str(debug_output)+" "\
                        +str(curr_start_seed)+" "\
                        +str(curr_end_seed)+" "\
                        +str(0)+" "\
                        +str(earliest_lineage_start_time)+" "\
                        +str(latest_lineage_start_time)+" "\
                        +str(events_end_time)+" "
        
        curr_list.append(validate_rate_command)
        command_count += 1
        
        for i in range(0,len(curr_list)):
            if deterministic_mode == 0:
                curr_list[i] = curr_list[i] + stochastic_theta_string
            if deterministic_mode == 1:
                curr_list[i] = curr_list[i] + deterministic_theta_string
                
        command_list = command_list + curr_list
        
    # Use processes equal to the number of cpus available
    cpu_count = multiprocessing.cpu_count()
    
    print(command_list)

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
