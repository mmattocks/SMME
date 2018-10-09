import multiprocessing
import os
import subprocess
import datetime
import numpy as np
from imageio.plugins._bsdf import BsdfSerializer

gomes_executable = '/home/main/chaste_build/projects/ISP/apps/GomesSimulator'
he_executable = '/home/main/chaste_build/projects/ISP/apps/HeSimulator'
boije_executable = '/home/main/chaste_build/projects/ISP/apps/BoijeSimulator'

empirical_data = '/home/main/git/chaste/projects/ISP/empirical_lineages.csv'

if not(os.path.isfile(gomes_executable)):
    raise Exception('Could not find executable: ' + gomes_executable)
if not(os.path.isfile(he_executable)):
    raise Exception('Could not find executable: ' + he_executable)
if not(os.path.isfile(boije_executable)):
    raise Exception('Could not find executable: ' + boije_executable)

#########################
# SIMULATION PARAMETERS
#########################

#Define start and end RNG seeds; determines:
#No. lineages per loss function run
#unique sequence of RNG results for each lineage

directory = "KolmogorovSequences"
output_mode = 2 #sequence sampler
debug_output = 0 #0=off;1=on
start_seed = 0
end_seed = 9999

##########################
#GLOBAL MODEL PARAMETERS
##########################

number_traversal_lineages = 60

############################
# SPECIFIC MODEL PARAMETERS 
############################

#GOMES MODEL
g_end_time = 480
g_normal_mean = 3.9716
g_normal_sigma = .32839
g_pPP = .055
g_pPD = .221
g_pBC = .128
g_pAC = .106
g_pMG = .028

g_string = gomes_executable + " " + directory + " Gomes " + str(output_mode) + " " + str(debug_output) + " " + str(start_seed) + " " + str(end_seed) + " " + str(g_end_time) + " " + str(g_normal_mean) + " " + str(g_normal_sigma) + " " + str(g_pPP) + " " + str(g_pPD) + " " + str(g_pBC) + " " + str(g_pAC) + " " + str(g_pMG)

#HE MODEL - GENERAL
h_fixture = 2
h_ath5founder = 0
h_start_time = 0
h_end_time = 80

#HE MODEL - STOCHASTIC
h_deterministic_mode = 0
h_mitotic_mode_phase_2 = 8 #These are phase lengths, so
h_mitotic_mode_phase_3 = 7 #Phase3 boundary = mmp2 + mmp3
h_phase_1_pPP = 1.0
h_phase_1_pPD = 0.0
h_phase_2_pPP = 0.2
h_phase_2_pPD = 0.4
h_phase_3_pPP = 0.2
h_phase_3_pPD = 0.0

h_string = he_executable + " " + directory + " He " + str(output_mode) + " " + str(h_deterministic_mode) + " " + str(h_fixture) + " " + str(h_ath5founder) + " " + str(debug_output) + " " + str(start_seed) + " " + str(end_seed) + " " + str(h_start_time) + " " + str(h_start_time) + " " + str(h_end_time) + " " + str(h_mitotic_mode_phase_2) + " " + str(h_mitotic_mode_phase_3) + " " +str(h_phase_1_pPP) + " " + str(h_phase_1_pPD) + " " + str(h_phase_2_pPP) + " " + str(h_phase_2_pPD) + " " + str(h_phase_3_pPP) + " " + str(h_phase_3_pPD)

#HE MODEL -STOCHASTIC (REFIT)
hr_mitotic_mode_phase_2 = 2.997 #These are phase lengths, so
hr_mitotic_mode_phase_3 = 9.730 #Phase3 boundary = mmp2 + mmp3
hr_phase_1_pPP = 1.0
hr_phase_1_pPD = 0.0
hr_phase_2_pPP = 0.2728
hr_phase_2_pPD = 0.4859
hr_phase_3_pPP = 0.2843
hr_phase_3_pPD = 0.0

hr_string = he_executable + " " + directory + " HeRefit " + str(output_mode) + " " + str(h_deterministic_mode) + " " + str(h_fixture) + " " + str(h_ath5founder) + " " + str(debug_output) + " " + str(start_seed) + " " + str(end_seed) + " " + str(h_start_time) + " " + str(h_start_time) + " " + str(h_end_time) + " " + str(hr_mitotic_mode_phase_2) + " " + str(hr_mitotic_mode_phase_3) + " " +str(hr_phase_1_pPP) + " " + str(hr_phase_1_pPD) + " " + str(hr_phase_2_pPP) + " " + str(hr_phase_2_pPD) + " " + str(hr_phase_3_pPP) + " " + str(hr_phase_3_pPD)

#HE MODEL - DETERMINISTIC ALTERNATIVE
d_deterministic_mode = 1
#Phase boundary shift parameters
d_phase_1_shape = 2.789
d_phase_1_scale = 1.430
d_phase_2_shape = 3.830
d_phase_2_scale = 3.524
d_phase_sister_shift_widths = .7288
d_phase_offset = .07611

d_string = he_executable + " " + directory + " Deterministic " + str(output_mode) + " " + str(d_deterministic_mode) + " " + str(h_fixture) + " " + str(h_ath5founder) + " " + str(debug_output) + " " + str(start_seed) + " " + str(end_seed) + " " + str(h_start_time) + " " + str(h_start_time) + " " + str(h_end_time) + " " + str(d_phase_1_shape) + " " + str(d_phase_1_scale) + " " + str(d_phase_2_shape) + " " + str(d_phase_2_scale) + " " + str(d_phase_sister_shift_widths) + " " + str(d_phase_offset)

#BOIJE MODEL
b_end_generation = 250
b_phase_2_generation = 3
b_phase_3_generation = 5
b_pAtoh7 = .32
b_pPtf1a = .3
b_png = .8

b_string = boije_executable + " " + directory + " Boije " + str(output_mode) + " " + str(debug_output) + " " + str(start_seed) + " " + str(end_seed) + " " + str(b_end_generation) + " " + str(b_phase_2_generation) + " " + str(b_phase_3_generation) + " " + str(b_pAtoh7) + " " + str(b_pPtf1a) + " " + str(b_png)


def main():

    command_list = []
    command_list.append(g_string)
    command_list.append(h_string)
    command_list.append(hr_string)
    command_list.append(d_string)
    command_list.append(b_string)
     
    # Use processes equal to the number of cpus available
    cpu_count = multiprocessing.cpu_count()
    
    print(command_list)

    print("Starting simulations with " + str(cpu_count) + " processes")
   
    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=cpu_count)

    # Pass the list of bash commands to the pool, block until pool is complete
    pool.map(execute_command, command_list, 1)
    
    e_sequences = traverse_lineages(empirical_data)
    g_sequences = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory + "/Gomes", skiprows=1, usecols=2)
    h_sequences = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory + "/He", skiprows=1, usecols=2)
    hr_sequences = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory + "/HeRefit", skiprows=1, usecols=2)
    d_sequences = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory + "/Deterministic", skiprows=1, usecols=2)
    b_sequences = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory + "/Boije", skiprows=1, usecols=2)

    e_acss = np.array(racss(e_sequences, alphabet=4))
    g_acss = np.array(racss(g_sequences, alphabet=4))
    h_acss = np.array(racss(h_sequences, alphabet=4))
    hr_acss = np.array(racss(hr_sequences, alphabet=4))
    d_acss = np.array(racss(d_sequences, alphabet=4))
    b_acss = np.array(racss(b_sequences, alphabet=4))

    e_mean_k = np.mean(e_acss[:,0])
    g_mean_k = np.mean(g_acss[:,0])
    h_mean_k = np.mean(h_acss[:,0])
    hr_mean_k = np.mean(hr_acss[:,0])
    d_mean_k = np.mean(d_acss[:,0])
    b_mean_k = np.mean(b_acss[:,0])


def traverse_lineages(filename):
    
    mode_sequences = np.zeros(end_seed+1)
    
    #load lineage tracing data
    e_data = np.loadtxt(filename, skiprows = 1, usecols = (0,1,2,3))

    p_RNG = np.random.RandomState(1)

    for curr_seed in range (start_seed,end_seed+1):
        k=0
        random_lineage_number = p_RNG.randint(1,number_traversal_lineages)
        lineage_events = np.array(e_data[np.where(e_data[:,0]==random_lineage_number)])
        
        #find mitotic mode of first event and write to log
        current_event = 1
        mitotic_mode = int(lineage_events[np.where(lineage_events[:,1]==current_event)][3])
        mode_sequences[k] = str(mitotic_mode)
        
        while mitotic_mode != 2:
            child_events = lineage_events[np.where(lineage_events[:,2]==current_event)]
            random_child_row = p_RNG.randint(0,np.ma.size(child_events,0))
            current_event= child_events[random_child_row,1]
            mitotic_mode = child_events[random_child_row,3]
            mode_sequences[k] += str(mitotic_mode)
            
        k += 1

    return mode_sequences

# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    print("Executing command: " + cmd)
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()



