import os
import numpy as np
import matplotlib.pyplot as plt
from fileinput import filename
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

count_seeds = 10000
event_seeds = 1000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations

#stats & line_plotter utility arrays
bin_sequence_24_32 = np.arange(1,26,1)
x_sequence_24_32 = np.arange(1,25,1)

bin_sequence_48 = np.arange(1,11,1)
x_sequence_48 = np.arange(1,10,1)

bin_sequence_wan = np.arange(1,16,1)
x_sequence_wan = np.arange(1,15,1)

bin_sequence_events = np.arange(30,85,5)
x_sequence_events = np.arange(30,80,5)

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

prob_empirical_24,bin_edges = np.histogram(np.sum(raw_counts_24,axis=1),bin_sequence_24_32,density=True) #obtain probability density histogram for counts, retabulating by summing across types
prob_empirical_32,bin_edges = np.histogram(np.sum(raw_counts_32,axis=1),bin_sequence_24_32,density=True)
prob_empirical_48,bin_edges = np.histogram(np.sum(raw_counts_48,axis=1),bin_sequence_48,density=True)

#no. of lineages observed per induction timepoint/event group
lineages_sampled_24 = 64
lineages_sampled_32 = 169
lineages_sampled_48 = 163
lineages_sampled_events = 60

#Mitotic mode rate probability arrays
raw_events = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_lineages.csv', skiprows=1, usecols=(3,5,8))
observed_events = raw_events[np.where(raw_events[:,2]==1)] #exclude any mitosis whose time was too early for recording

observed_PP = observed_events[np.where(observed_events[:,0]==0)]
observed_PD = observed_events[np.where(observed_events[:,0]==1)]
observed_DD = observed_events[np.where(observed_events[:,0]==2)]

prob_empirical_PP,bin_edges = np.histogram(observed_PP,bin_sequence_events,density=True)
prob_empirical_PD,bin_edges = np.histogram(observed_PD,bin_sequence_events,density=True)
prob_empirical_DD,bin_edges = np.histogram(observed_DD,bin_sequence_events,density=True)

#mutant results
wt_mean_clone_size = 6
ath5_mean_clone_size = 8
wt_even_odd_ratio = 1.1
ath5_even_odd_ratio = 2.2

##############################
# WAN ET AL EMPIRICAL RESULTS
##############################
#(derived from figure)
prob_wan = [0,0,.3143,.1857,.1757,.0871,.0871,.0786,.0271,.0179,.0071,0,0,.0071];
lineages_sampled_wan = 40 #not reported- approximation based on 95% CI

def main():
    #LOAD & PARSE HE MODEL OUTPUT FILES
    #original fit stochastic mitotic mode files
    o_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction24SDMode2", skiprows=1, usecols=3)
    o_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction32SDMode2", skiprows=1, usecols=3)
    o_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction48SDMode2", skiprows=1, usecols=3)
    o_events = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeModeEvent/inductionModeSDMode2", skiprows=1, usecols=(0,3))
    o_events_PP = np.array(o_events[np.where(o_events[:,1]==0)])
    o_events_PD = np.array(o_events[np.where(o_events[:,1]==1)])
    o_events_DD = np.array(o_events[np.where(o_events[:,1]==2)])
    o_counts_wan = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanSDMode2", skiprows=1, usecols=3)
    o_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanNoShadSDMode2", skiprows=1, usecols=3)
    o_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/ath5SDMode2", skiprows=1, usecols=3)
    o_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/validateSDMode2", skiprows=1, usecols=3)
    
    #refit stochastic mitotic mode files
    s_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction24SDMode0", skiprows=1, usecols=3)
    s_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction32SDMode0", skiprows=1, usecols=3)
    s_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/induction48SDMode0", skiprows=1, usecols=3)
    s_events = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeModeEvent/inductionModeSDMode0", skiprows=1, usecols=(0,3))
    s_events_PP = np.array(s_events[np.where(s_events[:,1]==0)])
    s_events_PD = np.array(s_events[np.where(s_events[:,1]==1)])
    s_events_DD = np.array(s_events[np.where(s_events[:,1]==2)])
    s_counts_wan = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanSDMode0", skiprows=1, usecols=3)
    s_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanNoShadSDMode0", skiprows=1, usecols=3)
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
    d_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/wanNoShadSDMode1", skiprows=1, usecols=3)
    d_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/ath5SDMode1", skiprows=1, usecols=3)
    d_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/HeCounts/validateSDMode1", skiprows=1, usecols=3)

    #calculate Ath5 clone size and even-odd rations
    o_wt_mean_clone_size = np.mean(o_counts_validate)
    o_ath5_mean_clone_size = np.mean(o_counts_ath5)

    s_wt_mean_clone_size = np.mean(s_counts_validate)
    s_ath5_mean_clone_size = np.mean(s_counts_ath5)
    
    d_wt_mean_clone_size = np.mean(d_counts_validate)
    d_ath5_mean_clone_size = np.mean(d_counts_ath5)
    
    o_wt_even_odd_ratio = len(np.where(o_counts_validate % 2 == 0)[0])/len(np.where(o_counts_validate % 2 == 1)[0])
    o_ath5_even_odd_ratio = len(np.where(o_counts_ath5 % 2 == 0)[0])/len(np.where(o_counts_ath5 % 2 == 1)[0])
    
    s_wt_even_odd_ratio = len(np.where(s_counts_validate % 2 == 0)[0])/len(np.where(s_counts_validate % 2 == 1)[0])
    s_ath5_even_odd_ratio = len(np.where(s_counts_ath5 % 2 == 0)[0])/len(np.where(s_counts_ath5 % 2 == 1)[0])

    d_wt_even_odd_ratio = len(np.where(d_counts_validate % 2 == 0)[0])/len(np.where(d_counts_validate % 2 == 1)[0])
    d_ath5_even_odd_ratio = len(np.where(d_counts_ath5 % 2 == 0)[0])/len(np.where(d_counts_ath5 % 2 == 1)[0])
    
    #SETUP FIGURES
    #4 figures to be generated; 3 show output of individual model modes, 4th overlays stochastic refit on deterministic fit
    o_fig, ((o_24,o_32,o_48),(o_PP,o_PD,o_DD),(o_wan, o_ath5size, o_ratio),(o_unused1,o_wan_noshad,o_unused2)) = plt.subplots(4,3,figsize=(12,8))
    s_fig, ((s_24,s_32,s_48),(s_PP,s_PD,s_DD),(s_wan, s_ath5size, s_ratio),(s_unused1,s_wan_noshad,s_unused2)) = plt.subplots(4,3,figsize=(12,8))
    d_fig, ((d_24,d_32,d_48),(d_PP,d_PD,d_DD),(d_wan, d_ath5size, d_ratio),(d_unused1,d_wan_noshad,d_unused2)) = plt.subplots(4,3,figsize=(12,8))
    sd_fig, ((sd_24,sd_32,sd_48),(sd_PP,sd_PD,sd_DD),(sd_wan, sd_ath5size, sd_ratio),(sd_unused1,sd_wan_noshad,sd_unused2)) = plt.subplots(4,3,figsize=(12,8))
    
    #turn off unused axes
    o_unused1.axis('off');o_unused2.axis('off');s_unused1.axis('off');s_unused2.axis('off'); d_unused1.axis('off'); d_unused2.axis('off'); sd_unused1.axis('off'); sd_unused2.axis('off');

    #set xlabels
    o_24.set_xlabel('Count');o_32.set_xlabel('Count');o_48.set_xlabel('Count');o_wan.set_xlabel('Count');o_wan_noshad.set_xlabel('Count')
    o_PP.set_xlabel('Age (hpf)');o_PD.set_xlabel('Age (hpf)');o_DD.set_xlabel('Age (hpf)')

    s_24.set_xlabel('Count');s_32.set_xlabel('Count');s_48.set_xlabel('Count');s_wan.set_xlabel('Count');s_wan_noshad.set_xlabel('Count')
    s_PP.set_xlabel('Age (hpf)');s_PD.set_xlabel('Age (hpf)');s_DD.set_xlabel('Age (hpf)')

    d_24.set_xlabel('Count');d_32.set_xlabel('Count');d_48.set_xlabel('Count');d_wan.set_xlabel('Count');d_wan_noshad.set_xlabel('Count')
    d_PP.set_xlabel('Age (hpf)');d_PD.set_xlabel('Age (hpf)');d_DD.set_xlabel('Age (hpf)')

    sd_24.set_xlabel('Count');sd_32.set_xlabel('Count');sd_48.set_xlabel('Count');sd_wan.set_xlabel('Count');sd_wan_noshad.set_xlabel('Count')
    sd_PP.set_xlabel('Age (hpf)');sd_PD.set_xlabel('Age (hpf)');sd_DD.set_xlabel('Age (hpf)')
    
    #set ylabels
    o_24.set_ylabel('Probability');o_32.set_ylabel('Probability');o_48.set_ylabel('Probability');o_wan.set_ylabel('Probability');o_wan_noshad.set_xlabel('Probability')
    o_PP.set_ylabel('PP Events/lineage/hr');o_PD.set_ylabel('PD Events/lineage/hr');o_DD.set_ylabel('DD Events/lineage/hr');
    o_ath5size.set_ylabel('Clone size');o_ratio.set_ylabel('Clone even/odd ratio')

    s_24.set_ylabel('Probability');s_32.set_ylabel('Probability');s_48.set_ylabel('Probability');s_wan.set_ylabel('Probability');s_wan_noshad.set_xlabel('Probability')
    s_PP.set_ylabel('PP Events/lineage/hr');s_PD.set_ylabel('PD Events/lineage/hr');s_DD.set_ylabel('DD Events/lineage/hr');
    s_ath5size.set_ylabel('Clone size');s_ratio.set_ylabel('Clone even/odd ratio')
    
    d_24.set_ylabel('Probability');d_32.set_ylabel('Probability');d_48.set_ylabel('Probability');d_wan.set_ylabel('Probability');d_wan_noshad.set_xlabel('Probability')
    d_PP.set_ylabel('PP Events/lineage/hr');d_PD.set_ylabel('PD Events/lineage/hr');d_DD.set_ylabel('DD Events/lineage/hr');
    d_ath5size.set_ylabel('Clone size');d_ratio.set_ylabel('Clone even/odd ratio')
    
    sd_24.set_ylabel('Probability');sd_32.set_ylabel('Probability');sd_48.set_ylabel('Probability');sd_wan.set_ylabel('Probability');sd_wan_noshad.set_xlabel('Probability')
    sd_PP.set_ylabel('PP Events/lineage/hr');sd_PD.set_ylabel('PD Events/lineage/hr');sd_DD.set_ylabel('DD Events/lineage/hr');
    sd_ath5size.set_ylabel('Clone size');sd_ratio.set_ylabel('Clone even/odd ratio')

    #PLOT DATA
    #plot lines, CIs, and bars
    line_plotter(o_24, o_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_32, o_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_48, o_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_PP, o_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_PD, o_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_DD, o_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_wan, o_counts_wan, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_wan_noshad, o_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'b-', '#000080', '#0000FF', 'y+', 0)
    o_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (o_wt_mean_clone_size/wt_mean_clone_size), (ath5_mean_clone_size/wt_mean_clone_size), (o_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#000080','#000080'],edgecolor='#000000',tick_label=['Observed WT','OSM WT', 'Observed Ath5mo', 'OSM Ath5mo'])
    o_ratio.bar([0,1,2,3], [wt_even_odd_ratio, o_wt_even_odd_ratio, ath5_even_odd_ratio, o_ath5_even_odd_ratio],color=['#000000','#000080','#000000','#000080'],edgecolor='#000000',tick_label=['Observed WT','OFit WT', 'Observed Ath5mo', 'OFit Ath5mo'])


    line_plotter(s_24, s_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_32, s_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_48, s_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_PP, s_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_PD, s_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_DD, s_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_wan, s_counts_wan, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_wan_noshad, s_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 0)
    s_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (s_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(s_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#800080', '#000000', '#800080'],edgecolor='#000000',tick_label=['Observed WT','SM WT', 'Observed Ath5mo', 'SM Ath5mo'])
    s_ratio.bar([0,1,2,3], [wt_even_odd_ratio, s_wt_even_odd_ratio, ath5_even_odd_ratio, s_ath5_even_odd_ratio],color=['#000000','#000080','#000000','#000080'],edgecolor='#000000',tick_label=['Observed WT','SM WT', 'Observed Ath5mo', 'SM Ath5mo'])
    
    line_plotter(d_24, d_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_32, d_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_48, d_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_PP, d_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_PD, d_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_DD, d_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_wan, d_counts_wan, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_wan_noshad, d_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 0)
    d_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (d_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(d_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#008000','#000000','#008000'],edgecolor='#000000',tick_label=['Observed WT','DM WT', 'Observed Ath5mo', 'DM Ath5mo'])
    d_ratio.bar([0,1,2,3], [wt_even_odd_ratio, d_wt_even_odd_ratio, ath5_even_odd_ratio, d_ath5_even_odd_ratio],color=['#000000','#008000','#000000','#008000'],edgecolor='#000000',tick_label=['Observed WT','DM WT', 'Observed Ath5mo', 'DM Ath5mo'])
    
    line_plotter(sd_24, s_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(sd_24, d_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(sd_32, s_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(sd_32, d_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(sd_48, s_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(sd_48, d_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(sd_PP, s_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(sd_PP, d_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(sd_PD, s_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(sd_PD, d_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(sd_DD, s_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(sd_DD, d_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(sd_wan, s_counts_wan, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(sd_wan, d_counts_wan, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(sd_wan_noshad, s_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(sd_wan_noshad, d_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 0)
    sd_ath5size.bar([0,1,2,3,4,5], [(wt_mean_clone_size/wt_mean_clone_size), (s_wt_mean_clone_size/wt_mean_clone_size),(d_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(s_ath5_mean_clone_size/wt_mean_clone_size),(d_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#800080','#008000','#000000','#800080','#008000'],edgecolor='#000000',tick_label=['Observed WT','SM WT','DM WT', 'Observed Ath5mo', 'SM Ath5Mo', 'DM Ath5mo'])
    sd_ratio.bar([0,1,2,3,4,5], [wt_even_odd_ratio, s_wt_even_odd_ratio, d_wt_even_odd_ratio, ath5_even_odd_ratio, s_ath5_even_odd_ratio, d_ath5_even_odd_ratio],color=['#000000','#800080','#008000','#000000','#800080','#008000'],edgecolor='#000000',tick_label=['Observed WT','SM WT','DM WT', 'Observed Ath5mo', 'SM Ath5Mo', 'DM Ath5mo'])

    #plot annotiations
    o_fig.patches.extend([plt.Rectangle((0.25,0.5),0.25,0.25,
                                  fill=True, color='g', alpha=0.5, zorder=1000,
                                  transform=o_fig.transFigure, figure=o_fig)])

    plt.tight_layout()
    plt.show()

#line plotter plots average counts or event event +- 2 standard deviations (He et al's "95% probability interval")
#"mode" 1 gives correct hourly output for 5-hour event bins
def line_plotter(subplot, data, bin_sequence, x_sequence, seeds, samples, empirical_prob, line_colour_string, fill_edge_colour_string, fill_face_colour_string, empirical_colour_string, mode):

    prob_histo, bin_edges = np.histogram(data, bin_sequence, density=True)
    
    #estimate of 95% CI for empirical observations of model output from repeated sampling of the data w/ the appropriate number of empirically observed lineages
    interval = sampler(data,samples,bin_sequence)
    subplot.plot(x_sequence,prob_histo, line_colour_string)
    subplot.fill_between(x_sequence, (prob_histo - interval), (prob_histo + interval), alpha=0.1, edgecolor=fill_edge_colour_string, facecolor=fill_face_colour_string)

    if mode == 0:
        subplot.plot(np.arange(1,len(empirical_prob)+1,1),empirical_prob, empirical_colour_string)
    
    if mode == 1:
        subplot.plot(np.arange(30,80,5),empirical_prob, empirical_colour_string)

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
    #these numpy arrays hold the individual RSS vals for timepoints + events
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

        
    events_plus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + "eventPlus", skiprows=1, usecols=(0,3))
    events_minus = np.loadtxt("/home/main/git/chaste/projects/ISP/testoutput/" + directory_name + "/" + file_name + "eventMinus", skiprows=1, usecols=(0,3))
    
    for i in range (0,3):
        mode_event_plus = np.array(events_plus[np.where(events_plus[:,1]==i)])
        mode_event_minus = np.array(events_minus[np.where(events_minus[:,1]==i)])

        histo_event_plus, bin_edges = np.histogram(mode_event_plus[:,0], event_bin_sequence, density=False)
        histo_event_minus, bin_edges = np.histogram(mode_event_minus[:,0], event_bin_sequence, density=False)
        
        prob_histo_event_plus = histo_event_plus / ((event_end_seed + 1)*5)
        prob_histo_event_minus = histo_event_minus / ((event_end_seed + 1)*5)

        residual_plus = prob_histo_event_plus - events_list[i]
        residual_minus = prob_histo_event_minus - events_list[i]
        
        rss_plus[i+3] = np.sum(np.square(residual_plus))
        rss_minus[i+3] = np.sum(np.square(residual_minus))
        
        plotter(plot_list[i+3], mode_event_plus[:,0], mode_event_minus[:,0], events_list[i],lineages_sampled_events, 1)

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
