import multiprocessing
import os
import subprocess
import datetime

import numpy as np
import matplotlib.pyplot as plt
from fileinput import filename
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

count_seeds = 100000
event_seeds = 10000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations


##############################
# HE ET AL EMPIRICAL RESULTS
##############################
#no. of lineages observed per induction timepoint/event group
lineages_sampled_24 = 64
lineages_sampled_32 = 169
lineages_sampled_48 = 163
lineages_sampled_events = 60

#binned histograms of counts at each timepoint
he_24_counts = np.array([ 0, 0, 1, 0, 1, 2, 0, 7, 4, 9, 6, 5, 5, 6, 5, 5, 1, 2, 2, 2, 0, 1 ])
he_32_counts = np.array([ 6, 20, 25, 22, 24, 17, 11, 7, 8, 5, 6, 3, 1, 6, 4, 1, 0, 0, 0, 2, 0, 0, 0, 1 ])
he_48_counts = np.array([ 59, 86, 2, 12, 1, 2, 0, 1 ])

prob_24_counts = np.array(he_24_counts/lineages_sampled_24)
prob_32_counts = np.array(he_32_counts/lineages_sampled_32)
prob_48_counts = np.array(he_48_counts/lineages_sampled_48)

#lineages giving mitotic mode events in 2 hr bins, starting at 30hr, from 60 empirically followed lineages
he_PP_events = np.array([8,23,16,10,3,1,0,0,0,0])
he_PD_events = np.array([2,21,7,10,2,0,0,0,0,0])
he_DD_events = np.array([0,3,13,26,33,29,16,3,4,0])

prob_empirical_PP = np.array((he_PP_events / lineages_sampled_events)/5)
prob_empirical_PD = np.array((he_PD_events / lineages_sampled_events)/5)
prob_empirical_DD = np.array((he_DD_events / lineages_sampled_events)/5)

#mutant results
wt_mean_clone_size = 6
ath5_mean_clone_size = 8
wt_even_odd_ratio = 1.1
ath5_even_odd_ratio = 2.2

#Wan et al. empirical results
prob_wan = [0,0,.3143,.1857,.1757,.0871,.0871,.0786,.0271,.0179,.0071,0,0,.0071];
lineages_sampled_wan = 40

def main():
    #load & parse He model output files
    #stochastic mitotic mode files
    s_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction24SDMode0", skiprows=1, usecols=3)
    s_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction32SDMode0", skiprows=1, usecols=3)
    s_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction48SDMode0", skiprows=1, usecols=3)
    s_events = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeModeEvent/inductionModeSDMode0", skiprows=1, usecols=(0,3))
    s_events_PP = np.array(s_events[np.where(s_events[:,1]==0)])
    s_events_PD = np.array(s_events[np.where(s_events[:,1]==1)])
    s_events_DD = np.array(s_events[np.where(s_events[:,1]==2)])
    s_counts_wan = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanSDMode0", skiprows=1, usecols=3)
    s_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/ath5SDMode0", skiprows=1, usecols=3)
    s_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/validateSDMode0", skiprows=1, usecols=3)    
    
    #deterministic mitotic mode files
    d_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction24SDMode1", skiprows=1, usecols=3)
    d_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction32SDMode1", skiprows=1, usecols=3)
    d_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction48SDMode1", skiprows=1, usecols=3)
    d_events = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeModeEvent/inductionModeSDMode1", skiprows=1, usecols=(0,3))
    d_events_PP = np.array(d_events[np.where(d_events[:,1]==0)])
    d_events_PD = np.array(d_events[np.where(d_events[:,1]==1)])
    d_events_DD = np.array(d_events[np.where(d_events[:,1]==2)])
    d_counts_wan = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanSDMode1", skiprows=1, usecols=3)
    d_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/ath5SDMode1", skiprows=1, usecols=3)
    d_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/validateSDMode1", skiprows=1, usecols=3)    

    s_wt_mean_clone_size = np.mean(s_counts_validate)
    s_ath5_mean_clone_size = np.mean(s_counts_ath5)
    
    d_wt_mean_clone_size = np.mean(d_counts_validate)
    d_ath5_mean_clone_size = np.mean(d_counts_ath5)
    
    s_wt_even_odd_ratio = len(np.where(s_counts_validate % 2 == 0)[0])/len(np.where(s_counts_validate % 2 == 1)[0])
    s_ath5_even_odd_ratio = len(np.where(s_counts_ath5 % 2 == 0)[0])/len(np.where(s_counts_ath5 % 2 == 1)[0])

    d_wt_even_odd_ratio = len(np.where(d_counts_validate % 2 == 0)[0])/len(np.where(d_counts_validate % 2 == 1)[0])
    d_ath5_even_odd_ratio = len(np.where(d_counts_ath5 % 2 == 0)[0])/len(np.where(d_counts_ath5 % 2 == 1)[0])
    
    #setup plot
    fig, ((ind24, ind32, ind48),(ratePP, ratePD, rateDD),(wan, ath5size, ath5ratio)) = plt.subplots(3,3,figsize=(12,6))
    
    #plot count lines
    line_plotter(ind24, s_counts_24, d_counts_24, np.arange(1,26,1), np.arange(1,25,1), count_seeds, lineages_sampled_24, prob_24_counts, 0)
    line_plotter(ind32, s_counts_32, d_counts_32, np.arange(1,26,1), np.arange(1,25,1), count_seeds, lineages_sampled_32, prob_32_counts, 0)
    line_plotter(ind48, s_counts_48, d_counts_48, np.arange(1,11,1), np.arange(1,10,1), count_seeds, lineages_sampled_48, prob_48_counts, 0)

    #plot events
    line_plotter(ratePP, s_events_PP[:,0], d_events_PP[:,0], np.arange(30,85,5), np.arange(30,80,5), event_seeds, lineages_sampled_events, prob_empirical_PP, 1)
    line_plotter(ratePD, s_events_PD[:,0], d_events_PD[:,0], np.arange(30,85,5), np.arange(30,80,5), event_seeds, lineages_sampled_events, prob_empirical_PD, 1)
    line_plotter(rateDD, s_events_DD[:,0], d_events_DD[:,0], np.arange(30,85,5), np.arange(30,80,5), event_seeds, lineages_sampled_events, prob_empirical_DD, 1)

    #plot Wan data
    line_plotter(wan, s_counts_wan, d_counts_wan, np.arange(1,16,1), np.arange(1,15,1), count_seeds, lineages_sampled_wan, prob_wan, 0)

    ath5size.bar([0,1,2,3,4,5],[(wt_mean_clone_size/wt_mean_clone_size),(s_wt_mean_clone_size/wt_mean_clone_size),(d_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(s_ath5_mean_clone_size/wt_mean_clone_size),(d_ath5_mean_clone_size/wt_mean_clone_size)])

    ath5ratio.bar([0,1,2,3,4,5],[wt_even_odd_ratio, s_wt_even_odd_ratio, d_wt_even_odd_ratio, ath5_even_odd_ratio, s_ath5_even_odd_ratio, d_ath5_even_odd_ratio])

    plt.show()

#data plotter function for monitoring SPSA results
def line_plotter(subplot,stochastic,deterministic, bin_sequence, x_sequence, seeds, samples, empirical_prob, mode):

    histo_s, bin_edges = np.histogram(stochastic, bin_sequence, density=False)
    histo_d, bin_edges = np.histogram(deterministic, bin_sequence, density=False)
    
    global prob_histo_s, prob_histo_d
    if mode == 0:
        prob_histo_s = histo_s / seeds
        prob_histo_d = histo_d / seeds
    
    if mode == 1:
        prob_histo_s = histo_s / (seeds*5)
        prob_histo_d = histo_d / (seeds*5)
    
    interval_s = sampler(stochastic,samples,bin_sequence)
    subplot.plot(x_sequence,prob_histo_s, 'm-')
    subplot.fill_between(x_sequence, (prob_histo_s - interval_s), (prob_histo_s + interval_s), alpha=0.1, edgecolor='#800080', facecolor='#FF00FF')

    interval_d = sampler(deterministic,samples,bin_sequence)
    subplot.plot(x_sequence,prob_histo_d, 'g-')
    subplot.fill_between(x_sequence, (prob_histo_d - interval_d), (prob_histo_d + interval_d), alpha=0.1, edgecolor='#008000', facecolor='#00FF00')

    if mode == 0:
        subplot.plot(np.arange(1,len(empirical_prob)+1,1),empirical_prob, 'k+')
    
    if mode == 1:
        subplot.plot(np.arange(30,80,5),empirical_prob, 'k+')

        
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



def evaluate_AIC_gradient(k, theta_plus, theta_minus, deterministic_mode, number_params, file_name):

    
    #these numpy arrays hold the individual RSS vals for timepoints + rates
    rss_plus = np.zeros(len(induction_times)+3)
    rss_minus = np.zeros(len(induction_times)+3)
    
    for i in range(0, len(plot_list)):
        plt.sca(plot_list[i])
        plt.cla()
    
    for i in range(0,len(induction_times)):
        counts_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Plus", skiprows=1, usecols=3)
        counts_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + str(induction_times[i]) + "Minus", skiprows=1, usecols=3)
        prob_histo_plus, bin_edges = np.histogram(counts_plus, bins=1000, range=(1,1001),density=True)
        prob_histo_minus, bin_edges = np.histogram(counts_minus, bins=1000, range=(1,1001),density=True)
        
        residual_plus = prob_histo_plus - empirical_prob_list[i]
        residual_minus = prob_histo_minus - empirical_prob_list[i]
        
        rss_plus[i] = np.sum(np.square(residual_plus))
        rss_minus[i] = np.sum(np.square(residual_minus))
        
        plotter(plot_list[i], counts_plus, counts_minus, empirical_prob_list[i],lineages_sampled_list[i], 0)

        
    rates_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + "RatePlus", skiprows=1, usecols=(0,3))
    rates_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + "RateMinus", skiprows=1, usecols=(0,3))
    
    for i in range (0,3):
        mode_rate_plus = np.array(rates_plus[np.where(rates_plus[:,1]==i)])
        mode_rate_minus = np.array(rates_minus[np.where(rates_minus[:,1]==i)])

        histo_rate_plus, bin_edges = np.histogram(mode_rate_plus[:,0], rate_bin_sequence, density=False)
        histo_rate_minus, bin_edges = np.histogram(mode_rate_minus[:,0], rate_bin_sequence, density=False)
        
        prob_histo_rate_plus = histo_rate_plus / ((rate_end_seed + 1)*5)
        prob_histo_rate_minus = histo_rate_minus / ((rate_end_seed + 1)*5)

        residual_plus = prob_histo_rate_plus - events_list[i]
        residual_minus = prob_histo_rate_minus - events_list[i]
        
        rss_plus[i+3] = np.sum(np.square(residual_plus))
        rss_minus[i+3] = np.sum(np.square(residual_minus))
        
        plotter(plot_list[i+3], mode_rate_plus[:,0], mode_rate_minus[:,0], events_list[i],lineages_sampled_events, 1)

    AIC_plus = 2 * number_params + number_comparisons * np.log(np.sum(rss_plus))
    AIC_minus = 2 * number_params + number_comparisons * np.log(np.sum(rss_minus))
    
    
    plt.show()
    
    log.write("theta_plus:\n")
    if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_plus[0]) + "\t" + str(theta_plus[1]) + "\t" + str(theta_plus[2]) + "\t" + str(theta_plus[3]) + "\t" + str(theta_plus[4]) + "\n")
    if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_plus[0]) + "\t" + str(theta_plus[1]) + "\t" + str(theta_plus[2]) + "\t" + str(theta_plus[3]) + "\t" + str(theta_plus[4]) + "\t" + str(theta_plus[5]) +"\n")
    log.write("theta_minus:\n")
    if deterministic_mode == 0: log.write(str(k) + "\t" + str(theta_minus[0]) + "\t" + str(theta_minus[1]) + "\t" + str(theta_minus[2]) + "\t" + str(theta_minus[3]) + "\t" + str(theta_minus[4]) + "\n")
    if deterministic_mode == 1: log.write(str(k) + "\t" + str(theta_minus[0]) + "\t" + str(theta_minus[1]) + "\t" + str(theta_minus[2]) + "\t" + str(theta_minus[3]) + "\t" + str(theta_minus[4]) + "\t" + str(theta_minus[5]) + "\n")

    log.write("PositiveAIC: " + str(AIC_plus) + " NegativeAIC: " + str(AIC_minus) + "\n")

    AIC_gradient_sample = AIC_plus - AIC_minus;

    return AIC_gradient_sample

# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
