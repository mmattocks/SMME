import os
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from io import BytesIO

#PLoS formatting stuff
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

#AIC & Plotting utility params
event_seeds = 1000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations

bin_sequence_events = np.arange(30,85,5)
x_sequence_events = np.arange(30,80,5)

##############################
# HE ET AL EMPIRICAL RESULTS
##############################

rate_bin_sequence = np.arange(30,85,5)
rate_x_sequence = np.arange(30,80,5)
rate_trim_value = 10

#no. of lineages EO per induction timepoint/event group
lineages_sampled_events = 60

#Mitotic mode rate probability arrays
raw_events = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_lineages.csv', skiprows=1, usecols=(3,5,8))
EO_events = raw_events[np.where(raw_events[:,2]==1)] #exclude any mitosis whose time was too early for recording

event_counts, bin_edges = np.histogram(EO_events,bin_sequence_events,density=False)
#hourly values
empirical_events_per_lineage = np.array(event_counts/(lineages_sampled_events*5))

def main():
    #LOAD & PARSE HE MODEL OUTPUT FILES
    #original fit stochastic mitotic mode files
    o_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode2", skiprows=1, usecols=(0))
    
    #refit stochastic mitotic mode files
    s_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode0", skiprows=1, usecols=(0))
     
    #deterministic mitotic mode files
    d_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode1", skiprows=1, usecols=(0))

    original = line_plotter(o_events,bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events, 'b--', '#000080', '#0000FF')
    refit = line_plotter(s_events,bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events, 'm-', '#800080', '#FF00FF')
    deterministic = line_plotter(d_events,bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events, 'g-', '#008000', '#00FF00')
    
    empirical = plt.plot(np.arange(30,80,5),empirical_events_per_lineage, 'k+')
    
    plt.legend((original[0], refit[0], deterministic[0], empirical[0]), ('SM (He fit)', 'SM (SPSA fit)', 'DM', 'Observed'))
    
    plt.xlabel ('Age (hpf)')
    plt.ylabel ('Probability of mitotic event per lineage per hour')
    
    
    plt.savefig("per_lineage_mitotic_rate.png")
    png_memory = BytesIO()
    plt.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('per_lineage_mitotic_rate.tiff')
    png_memory.close()
    
    plt.show()

#line plotter plots average counts or event event +- 2 standard deviations (He et al's "95% probability interval")
#"mode" 1 gives correct hourly output for 5-hour event bins. mode 2 begins Wan plots at clone size 2, reflecting unavailable data
def line_plotter(data, bin_sequence, x_sequence, seeds, samples, line_colour_string, fill_edge_colour_string, fill_face_colour_string):

    histo, bin_edges = np.histogram(data, bin_sequence, density=False)
    prob_histo = np.array(histo / (seeds * 5))
    
    #estimate of 95% CI for empirical observations of model output from repeated sampling of the data w/ the appropriate number of empirically observed lineages
    interval = sampler(data,samples,bin_sequence)
    line = plt.plot(x_sequence,prob_histo, line_colour_string)
    plt.fill_between(x_sequence, (prob_histo - interval), (prob_histo + interval), alpha=0.1, edgecolor=fill_edge_colour_string, facecolor=fill_face_colour_string)

    return line

def sampler(data,samples,bin_sequence):
    
    if len(data) == 0: #catches edge case w/ no entries for a mitotic mode
        data=[0]
    
    base_sample=np.zeros((error_samples,len(bin_sequence)-1))

    for i in range(0,error_samples):
        new_data_sample = np.random.choice(data,samples)
        new_histo, bin_edges = np.histogram(new_data_sample, bin_sequence, density=False)
        new_histo_prob = np.array(new_histo/samples)
        base_sample[i,:] = new_histo_prob
        
    sample_95CI = np.array(2 * (np.std(base_sample,0)))
    
    return sample_95CI

if __name__ == "__main__":
    main()
