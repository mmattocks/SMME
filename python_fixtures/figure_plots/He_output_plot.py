import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fileinput import filename
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

from PIL import Image
from io import BytesIO

#PLoS formatting stuff
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

#AIC & Plotting utility params
count_seeds = 10000
event_seeds = 1000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations

he_model_params = 15
det_model_params = 13

#stats & line_plotter utility arrays
bin_sequence_24_32 = np.arange(1,26,1)
x_sequence_24_32 = np.arange(1,25,1)

bin_sequence_48 = np.arange(1,11,1)
x_sequence_48 = np.arange(1,10,1)

bin_sequence_wan = np.arange(2,17,1)
x_sequence_wan = np.arange(2,16,1)

bin_sequence_events = np.arange(30,85,5)
x_sequence_events = np.arange(30,80,5)

##############################
# HE ET AL EMPIRICAL RESULTS
##############################

#Count probability arrays
raw_counts = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_counts.csv', skiprows=1, usecols=(3,4,5,6,7,8,9,10)) #collect the per-cell-type counts
raw_counts_24 = raw_counts[0:64,:]
raw_counts_32 = raw_counts[64:233,:]
raw_counts_48 = raw_counts[233:396,:]

prob_empirical_24,bin_edges = np.histogram(np.sum(raw_counts_24,axis=1),bin_sequence_24_32,density=True) #obtain probability density histogram for counts, retabulating by summing across types
prob_empirical_32,bin_edges = np.histogram(np.sum(raw_counts_32,axis=1),bin_sequence_24_32,density=True)
prob_empirical_48,bin_edges = np.histogram(np.sum(raw_counts_48,axis=1),bin_sequence_48,density=True)

#no. of lineages EO per induction timepoint/event group
lineages_sampled_24 = 64
lineages_sampled_32 = 169
lineages_sampled_48 = 163
lineages_sampled_events = 60

#Mitotic mode rate probability arrays
raw_events = np.loadtxt('/home/main/git/chaste/projects/ISP/empirical_data/empirical_lineages.csv', skiprows=1, usecols=(3,5,8))
EO_events = raw_events[np.where(raw_events[:,2]==1)] #exclude any mitosis whose time was too early for recording

EO_PP = EO_events[np.where(EO_events[:,0]==0)]
EO_PD = EO_events[np.where(EO_events[:,0]==1)]
EO_DD = EO_events[np.where(EO_events[:,0]==2)]

prob_empirical_PP,bin_edges = np.histogram(EO_PP,bin_sequence_events,density=True)
prob_empirical_PD,bin_edges = np.histogram(EO_PD,bin_sequence_events,density=True)
prob_empirical_DD,bin_edges = np.histogram(EO_DD,bin_sequence_events,density=True)

empirical_fit_list = [prob_empirical_24,prob_empirical_32,prob_empirical_48,prob_empirical_PP,prob_empirical_PD,prob_empirical_DD]

#mutant results
wt_mean_clone_size = 6
ath5_mean_clone_size = 8
wt_even_odd_ratio = 1.1
ath5_even_odd_ratio = 2.2

##############################
# WAN ET AL EMPIRICAL RESULTS
##############################
#(derived from figure)
prob_wan = [.3143,.1857,.1757,.0871,.0871,.0786,.0271,.0179,.0071,0,0,.0071];
lineages_sampled_wan = 40 #not reported- approximation based on 95% CI

empirical_test_list = [prob_wan, wt_mean_clone_size, ath5_mean_clone_size, wt_even_odd_ratio, ath5_even_odd_ratio]

def main():
    #LOAD & PARSE HE MODEL OUTPUT FILES
    #original fit stochastic mitotic mode files
    o_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction24SDMode2", skiprows=1, usecols=3)
    o_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction32SDMode2", skiprows=1, usecols=3)
    o_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction48SDMode2", skiprows=1, usecols=3)
    o_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode2", skiprows=1, usecols=(0,3))
    o_events_PP = np.array(o_events[np.where(o_events[:,1]==0)])
    o_events_PD = np.array(o_events[np.where(o_events[:,1]==1)])
    o_events_DD = np.array(o_events[np.where(o_events[:,1]==2)])
    o_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/wanNoShadSDMode2", skiprows=1, usecols=3)
    o_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/ath5SDMode2", skiprows=1, usecols=3)
    o_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/validateSDMode2", skiprows=1, usecols=3)
    
    #refit stochastic mitotic mode files
    s_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction24SDMode0", skiprows=1, usecols=3)
    s_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction32SDMode0", skiprows=1, usecols=3)
    s_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction48SDMode0", skiprows=1, usecols=3)
    s_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode0", skiprows=1, usecols=(0,3))
    s_events_PP = np.array(s_events[np.where(s_events[:,1]==0)])
    s_events_PD = np.array(s_events[np.where(s_events[:,1]==1)])
    s_events_DD = np.array(s_events[np.where(s_events[:,1]==2)])
    s_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/wanNoShadSDMode0", skiprows=1, usecols=3)
    s_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/ath5SDMode0", skiprows=1, usecols=3)
    s_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/validateSDMode0", skiprows=1, usecols=3)
    
    #deterministic mitotic mode files
    d_counts_24 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction24SDMode1", skiprows=1, usecols=3)
    d_counts_32 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction32SDMode1", skiprows=1, usecols=3)
    d_counts_48 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/induction48SDMode1", skiprows=1, usecols=3)
    d_events = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeModeEvent/inductionModeSDMode1", skiprows=1, usecols=(0,3))
    d_events_PP = np.array(d_events[np.where(d_events[:,1]==0)])
    d_events_PD = np.array(d_events[np.where(d_events[:,1]==1)])
    d_events_DD = np.array(d_events[np.where(d_events[:,1]==2)])
    d_counts_wan_no_shadow = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/wanNoShadSDMode1", skiprows=1, usecols=3)
    d_counts_ath5 = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/ath5SDMode1", skiprows=1, usecols=3)
    d_counts_validate = np.loadtxt("/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/HeCounts/validateSDMode1", skiprows=1, usecols=3)

    #calculate Ath5 clone size and even-odd ratios
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
    
    o_fit_output_list = [o_counts_24, o_counts_32, o_counts_48, o_events_PP, o_events_PD, o_events_DD]
    o_test_output_list = [o_counts_wan_no_shadow, o_wt_mean_clone_size, o_ath5_mean_clone_size, o_wt_even_odd_ratio, o_ath5_even_odd_ratio]
    
    s_fit_output_list = [s_counts_24, s_counts_32, s_counts_48, s_events_PP, s_events_PD, s_events_DD]
    s_test_output_list = [s_counts_wan_no_shadow, s_wt_mean_clone_size, s_ath5_mean_clone_size, s_wt_even_odd_ratio, s_ath5_even_odd_ratio]
    
    d_fit_output_list = [d_counts_24, d_counts_32, d_counts_48, d_events_PP, d_events_PD, d_events_DD]
    d_test_output_list = [d_counts_wan_no_shadow, d_wt_mean_clone_size, d_ath5_mean_clone_size, d_wt_even_odd_ratio, d_ath5_even_odd_ratio]
    
    #SETUP FIGURES
    #4 figures to be generated; 3 show output of individual model modes, 4th overlays stochastic refit on deterministic fit
    o_fig, ((o_24,o_32),(o_48,o_PP),(o_PD,o_DD),(o_wan_noshad, o_ath5size),(o_ratio,o_unused)) = plt.subplots(5,2,figsize=(7.5,8.75))
    s_fig, ((s_24,s_32),(s_48,s_PP),(s_PD,s_DD),(s_wan_noshad, s_ath5size),(s_ratio,s_unused)) = plt.subplots(5,2,figsize=(7.5,8.75))
    d_fig, ((d_24,d_32),(d_48,d_PP),(d_PD,d_DD),(d_wan_noshad, d_ath5size),(d_ratio,d_unused)) = plt.subplots(5,2,figsize=(7.5,8.75))
    sd_fig, ((sd_24,sd_32),(sd_48,sd_PP),(sd_PD,sd_DD),(sd_wan_noshad, sd_ath5size),(sd_ratio,sd_unused)) = plt.subplots(5,2,figsize=(7.5,8.75))

    #turn off unused axes
    o_unused.axis('off');s_unused.axis('off');d_unused.axis('off');sd_unused.axis('off');

    #set xlabels
    o_24.set_xlabel('Count');o_32.set_xlabel('Count');o_48.set_xlabel('Count');o_wan_noshad.set_xlabel('Count')
    o_PP.set_xlabel('Age (hpf)');o_PD.set_xlabel('Age (hpf)');o_DD.set_xlabel('Age (hpf)')

    s_24.set_xlabel('Count');s_32.set_xlabel('Count');s_48.set_xlabel('Count');s_wan_noshad.set_xlabel('Count')
    s_PP.set_xlabel('Age (hpf)');s_PD.set_xlabel('Age (hpf)');s_DD.set_xlabel('Age (hpf)')

    d_24.set_xlabel('Count');d_32.set_xlabel('Count');d_48.set_xlabel('Count');d_wan_noshad.set_xlabel('Count')
    d_PP.set_xlabel('Age (hpf)');d_PD.set_xlabel('Age (hpf)');d_DD.set_xlabel('Age (hpf)')

    sd_24.set_xlabel('Count');sd_32.set_xlabel('Count');sd_48.set_xlabel('Count');sd_wan_noshad.set_xlabel('Count')
    sd_PP.set_xlabel('Age (hpf)');sd_PD.set_xlabel('Age (hpf)');sd_DD.set_xlabel('Age (hpf)')
     
    #set ylabels
    o_24.set_ylabel('Probability');o_32.set_ylabel('Probability');o_48.set_ylabel('Probability');o_wan_noshad.set_ylabel('Probability')
    o_PP.set_ylabel('PP Probability');o_PD.set_ylabel('PD Probability');o_DD.set_ylabel('DD Probability');
    o_ath5size.set_ylabel('Clone size');o_ratio.set_ylabel('Clone even/odd ratio')
  
    s_24.set_ylabel('Probability');s_32.set_ylabel('Probability');s_48.set_ylabel('Probability');s_wan_noshad.set_ylabel('Probability')
    s_PP.set_ylabel('PP Probability');s_PD.set_ylabel('PD Probability');s_DD.set_ylabel('DD Probability');
    s_ath5size.set_ylabel('Clone size');s_ratio.set_ylabel('Clone even/odd ratio')
      
    d_24.set_ylabel('Probability');d_32.set_ylabel('Probability');d_48.set_ylabel('Probability');d_wan_noshad.set_ylabel('Probability')
    d_PP.set_ylabel('PP Probability');d_PD.set_ylabel('PD Probability');d_DD.set_ylabel('DD Probability');
    d_ath5size.set_ylabel('Clone size');d_ratio.set_ylabel('Clone even/odd ratio')
     
    sd_24.set_ylabel('Probability');sd_32.set_ylabel('Probability');sd_48.set_ylabel('Probability');sd_wan_noshad.set_ylabel('Probability')
    sd_PP.set_ylabel('PP Probability');sd_PD.set_ylabel('PD Probability');sd_DD.set_ylabel('DD Probability');
    sd_ath5size.set_ylabel('Clone size');sd_ratio.set_ylabel('Clone even/odd ratio')

    #PLOT DATA
    #plot lines, CIs, and bars
    o_objects = line_plotter(o_24, o_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_32, o_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_48, o_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'b-', '#000080', '#0000FF', 'y+', 0)
    line_plotter(o_PP, o_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_PD, o_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_DD, o_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'b-', '#000080', '#0000FF', 'y+', 1)
    line_plotter(o_wan_noshad, o_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'b-', '#000080', '#0000FF', 'y+', 2)
    o_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (o_wt_mean_clone_size/wt_mean_clone_size), (ath5_mean_clone_size/wt_mean_clone_size), (o_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#000080','#000080'],edgecolor='#000000',tick_label=['EO WT','OSM WT', 'EO Ath5mo', 'OSM Ath5mo'])
    o_ratio.bar([0,1,2,3], [wt_even_odd_ratio, o_wt_even_odd_ratio, ath5_even_odd_ratio, o_ath5_even_odd_ratio],color=['#000000','#000080','#000000','#000080'],edgecolor='#000000',tick_label=['EO WT','OFit WT', 'EO Ath5mo', 'OFit Ath5mo'])
    plt.setp(o_ath5size.get_xticklabels(), rotation=25, ha='right')
    plt.setp(o_ratio.get_xticklabels(), rotation=25, ha='right')

    s_objects = line_plotter(s_24, s_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_32, s_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_48, s_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'm-', '#800080', '#FF00FF', 'k+', 0)
    line_plotter(s_PP, s_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_PD, s_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_DD, s_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'m-', '#800080', '#FF00FF', 'k+', 1)
    line_plotter(s_wan_noshad, s_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 2)
    s_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (s_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(s_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#800080', '#000000', '#800080'],edgecolor='#000000',tick_label=['EO WT','SM WT', 'EO Ath5mo', 'SM Ath5mo'])
    s_ratio.bar([0,1,2,3], [wt_even_odd_ratio, s_wt_even_odd_ratio, ath5_even_odd_ratio, s_ath5_even_odd_ratio],color=['#000000','#000080','#000000','#000080'],edgecolor='#000000',tick_label=['EO WT','SM WT', 'EO Ath5mo', 'SM Ath5mo'])
    plt.setp(s_ath5size.get_xticklabels(), rotation=25, ha='right')
    plt.setp(s_ratio.get_xticklabels(), rotation=25, ha='right')

    
    d_objects = line_plotter(d_24, d_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_32, d_counts_32, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_32, prob_empirical_32, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_48, d_counts_48, bin_sequence_48, x_sequence_48, count_seeds, lineages_sampled_48, prob_empirical_48, 'g-', '#008000', '#00FF00', 'k+', 0)
    line_plotter(d_PP, d_events_PP[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PP,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_PD, d_events_PD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_PD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_DD, d_events_DD[:,0], bin_sequence_events, x_sequence_events, event_seeds, lineages_sampled_events,prob_empirical_DD,'g-', '#008000', '#00FF00', 'k+', 1)
    line_plotter(d_wan_noshad, d_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 2)
    d_ath5size.bar([0,1,2,3], [(wt_mean_clone_size/wt_mean_clone_size), (d_wt_mean_clone_size/wt_mean_clone_size),(ath5_mean_clone_size/wt_mean_clone_size),(d_ath5_mean_clone_size/wt_mean_clone_size)],color=['#000000','#008000','#000000','#008000'],edgecolor='#000000',tick_label=['EO WT','DM WT', 'EO Ath5mo', 'DM Ath5mo'])
    d_ratio.bar([0,1,2,3], [wt_even_odd_ratio, d_wt_even_odd_ratio, ath5_even_odd_ratio, d_ath5_even_odd_ratio],color=['#000000','#008000','#000000','#008000'],edgecolor='#000000',tick_label=['EO WT','DM WT', 'EO Ath5mo', 'DM Ath5mo'])
    plt.setp(d_ath5size.get_xticklabels(), rotation=25, ha='right')
    plt.setp(d_ratio.get_xticklabels(), rotation=25, ha='right')

    
    sd_objects_s = line_plotter(sd_24, s_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'm-', '#800080', '#FF00FF', 'k+', 0)
    sd_objects_d = line_plotter(sd_24, d_counts_24, bin_sequence_24_32, x_sequence_24_32, count_seeds, lineages_sampled_24, prob_empirical_24, 'g-', '#008000', '#00FF00', 'k+', 0)
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
    line_plotter(sd_wan_noshad, s_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'm-', '#800080', '#FF00FF', 'k+', 2)
    line_plotter(sd_wan_noshad, d_counts_wan_no_shadow, bin_sequence_wan, x_sequence_wan, count_seeds, lineages_sampled_wan, prob_wan, 'g-', '#008000', '#00FF00', 'k+', 2)
    sd_ath5size.bar([0,1,2,3,4,5], [(wt_mean_clone_size), (s_wt_mean_clone_size),(d_wt_mean_clone_size),(ath5_mean_clone_size),(s_ath5_mean_clone_size),(d_ath5_mean_clone_size)],color=['#000000','#800080','#008000','#000000','#800080','#008000'],edgecolor='#000000',tick_label=['EO WT','SM WT','DM WT', 'EO Ath5mo', 'SM Ath5Mo', 'DM Ath5mo'])
    sd_ratio.bar([0,1,2,3,4,5], [wt_even_odd_ratio, s_wt_even_odd_ratio, d_wt_even_odd_ratio, ath5_even_odd_ratio, s_ath5_even_odd_ratio, d_ath5_even_odd_ratio],color=['#000000','#800080','#008000','#000000','#800080','#008000'],edgecolor='#000000',tick_label=['EO WT','SM WT','DM WT', 'EO Ath5mo', 'SM Ath5Mo', 'DM Ath5mo'])
    plt.setp(sd_ath5size.get_xticklabels(), rotation=25, ha='right')
    plt.setp(sd_ratio.get_xticklabels(), rotation=25, ha='right')

    #calculate AICs, write table to file
    o_fit_AIC = calculate_fit_AIC(o_fit_output_list, empirical_fit_list, he_model_params)
    o_test_AIC = calculate_test_AIC(o_test_output_list, empirical_fit_list, he_model_params)
    
    s_fit_AIC = calculate_fit_AIC(s_fit_output_list, empirical_fit_list, he_model_params)
    s_test_AIC = calculate_test_AIC(s_test_output_list, empirical_test_list, he_model_params)
    
    d_fit_AIC = calculate_fit_AIC(d_fit_output_list, empirical_fit_list, det_model_params)
    d_test_AIC = calculate_test_AIC(d_test_output_list, empirical_test_list, det_model_params)
    
    AIC_table('AIC_table.tex', o_fit_AIC, o_test_AIC, s_fit_AIC, s_test_AIC, d_fit_AIC, d_test_AIC)
    
    #plot annotiations
    plot_annotater(s_fig,s_objects,'Stochastic Model')
    plot_annotater(o_fig,o_objects,'Stochastic Model (original fit)')
    plot_annotater(d_fig,d_objects, 'Deterministic Model')
    plot_annotater(sd_fig,sd_objects_s, 'Stochastic Model', sd_objects_d, 'Deterministic Model')
    
    #export as .tiff, .png
    o_fig.savefig('o_fig.png',dpi=600)
    png_memory = BytesIO()
    o_fig.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('o_fig.tiff')
    png_memory.close()
    
    s_fig.savefig('s_fig.png',dpi=600)
    png_memory = BytesIO()
    s_fig.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('s_fig.tiff')
    png_memory.close()
    
    d_fig.savefig('d_fig.png',dpi=600)
    png_memory = BytesIO()
    d_fig.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('d_fig.tiff')
    png_memory.close()

    sd_fig.savefig('sd_fig.png',dpi=600)
    png_memory = BytesIO()
    sd_fig.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('sd_fig.tiff')
    png_memory.close()
    
    plt.show()


def calculate_fit_AIC(output_list, empirical_list, number_params):
    #function calculates AIC for He et al counts & event probability density (NOT per-lineage event probability)
    #this numpy array holds the individual RSS vals for timepoints + events
    rss = np.zeros(6)
    number_comparisons = 0
        
    for i in range(0,3):
        output_counts = output_list[i]
        empirical_prob = empirical_list[i]
        
        output_histo,bin_edges = np.histogram(output_counts, np.arange(1,len(empirical_prob)+2), density = False)
        output_prob = np.array(output_histo/count_seeds)
        
        residual = np.array(output_prob- empirical_prob)
        number_comparisons += len(residual)
        rss[i] = np.sum(np.square(residual))
    
    for i in range (0,3):
        output_events = output_list[i+3]
        empirical_prob = empirical_list[i+3]
        
        histo_events,bin_edges = np.histogram(output_events, bin_sequence_events, density=True)
        output_prob = np.array(histo_events/event_seeds)

        residual = np.array(output_prob - empirical_prob)
        number_comparisons += len(residual)
        rss[i+3] = np.sum(np.square(residual))
        
    AIC = 2 * number_params + number_comparisons * np.log(np.sum(rss))
    
    return AIC

def calculate_test_AIC(output_list, empirical_list, number_params):
    #function calculates AIC for He mutant data & Wan counts
    rss=np.zeros(3)
    number_comparisons = 0
    
    #Wan data
    for i in range (0,1):
        output_counts = output_list[i]
        empirical_prob = empirical_list[i]
        
        output_histo,bin_edges = np.histogram(output_counts, np.arange(2,len(empirical_prob)+3), density = False)
        output_prob = np.array(output_histo/count_seeds)
        
        residual = np.array(output_prob- empirical_prob)
        number_comparisons += len(residual)
        rss[i] = np.sum(np.square(residual))
    
    #He data
    for i in range(1,3):
        output_measure = output_list[i]
        empirical_measure = empirical_list[i]
        
        residual = np.array(output_measure - empirical_measure)
        number_comparisons += 1
        rss[i] = np.sum(np.square(residual))
        
    AIC = 2 * number_params + number_comparisons * np.log(np.sum(rss))
    
    return AIC

#line plotter plots average counts or event event +- 2 standard deviations (He et al's "95% probability interval")
#"mode" 1 gives correct hourly output for 5-hour event bins. mode 2 begins Wan plots at clone size 2, reflecting unavailable data
def line_plotter(subplot, data, bin_sequence, x_sequence, seeds, samples, empirical_prob, line_colour_string, fill_edge_colour_string, fill_face_colour_string, empirical_colour_string, mode):

    legend_objects = []
    
    prob_histo, bin_edges = np.histogram(data, bin_sequence, density=True)
    
    #estimate of 95% CI for empirical observations of model output from repeated sampling of the data w/ the appropriate number of empirically observed lineages
    interval = sampler(data,samples,bin_sequence)
    legend_objects.append(subplot.plot(x_sequence,prob_histo, line_colour_string))
    subplot.fill_between(x_sequence, (prob_histo - interval), (prob_histo + interval), alpha=0.1, edgecolor=fill_edge_colour_string, facecolor=fill_face_colour_string)

    if mode == 0:
        legend_objects.append(subplot.plot(np.arange(1,len(empirical_prob)+1,1),empirical_prob, empirical_colour_string))
    
    if mode == 1:
        legend_objects.append(subplot.plot(np.arange(30,80,5),empirical_prob, empirical_colour_string))
        
    if mode == 2:
        legend_objects.append(subplot.plot(np.arange(2,len(empirical_prob)+2,1), empirical_prob, empirical_colour_string))
        
    return legend_objects

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

def plot_annotater(fig_to_ann, legend_objects_1, legend_objects_1_name, legend_objects_2=None, legend_objects_2_name=None):
    #subplot labels
    labels = ['A','B','C','D','E','F','G','H','I','']
    
    for i, ax in enumerate(fig_to_ann.axes):
        ax.text(.01,.95, labels[i], transform=ax.transAxes, fontweight='bold', va='top')

        time_rectangle_x = .70
        time_rectangle_y = .65
        time_rectangle_width = .2
        time_rectangle_height = .03
        time_caret_half_width = .02
        time_caret_height = np.sqrt(3)*(time_caret_half_width*2)
        
        if i < 3:
            #add the time bar and label 72 hr caret on appropriate plots
            ax.add_patch(plt.Rectangle((time_rectangle_x,time_rectangle_y),time_rectangle_width,time_rectangle_height, color='k', transform=ax.transAxes))
            ax.add_patch(plt.Polygon(((time_rectangle_x+time_rectangle_width,time_rectangle_y+time_rectangle_height),(time_rectangle_x+time_rectangle_width-time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height),(time_rectangle_x+time_rectangle_width+time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height)),color='k', fill=False, transform=ax.transAxes))
            ax.text(time_rectangle_x+time_rectangle_width,time_rectangle_y+time_rectangle_height+time_caret_height,'72hpf', transform=ax.transAxes, fontsize='8',va='bottom',ha='center')
        
        if i == 0:
            #add 24hr caret
            ax.add_patch(plt.Polygon(((time_rectangle_x,time_rectangle_y+time_rectangle_height),(time_rectangle_x-time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height),(time_rectangle_x+time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height)),color='k', fill=False, transform=ax.transAxes))
            ax.text(time_rectangle_x,time_rectangle_y+time_rectangle_height+time_caret_height,'24hpf', transform=ax.transAxes, fontsize='8',va='bottom',ha='center')
        
        if i == 1:
            #add 32hr caret
            ax.add_patch(plt.Polygon(((time_rectangle_x+(8/48)*time_rectangle_width,time_rectangle_y+time_rectangle_height),(time_rectangle_x+(8/48)*time_rectangle_width-time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height),(time_rectangle_x+(8/48)*time_rectangle_width+time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height)),color='k', fill=False, transform=ax.transAxes))
            ax.text(time_rectangle_x+(8/48)*time_rectangle_width,time_rectangle_y+time_rectangle_height+time_caret_height,'32hpf', transform=ax.transAxes, fontsize='8',va='bottom',ha='center')

        if i == 2:
            ax.add_patch(plt.Polygon(((time_rectangle_x+(24/48)*time_rectangle_width,time_rectangle_y+time_rectangle_height),(time_rectangle_x+(24/48)*time_rectangle_width-time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height),(time_rectangle_x+(24/48)*time_rectangle_width+time_caret_half_width,time_rectangle_y+time_rectangle_height+time_caret_height)),color='k', fill=False, transform=ax.transAxes))
            ax.text(time_rectangle_x+(24/48)*time_rectangle_width,time_rectangle_y+time_rectangle_height+time_caret_height,'48hpf', transform=ax.transAxes, fontsize='8',va='bottom',ha='center')
        
        parent_circle_coord = (.85,.8)
        child_1_coord = (.8,.6)
        child_2_coord = (.9,.6)
        ellipse_width = .03
        ellipse_height = 2.5*ellipse_width
        
        if i == 3:
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_1_coord),color='k',transform=ax.transAxes))
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_2_coord),color='k',transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(parent_circle_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_1_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_2_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.text(.91, .7, 'PP', transform=ax.transAxes, va='center', ha='left')

        if i == 4:
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_1_coord),color='k',transform=ax.transAxes))
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_2_coord),color='k',transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(parent_circle_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_1_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_2_coord,ellipse_width,ellipse_height, edgecolor='k', facecolor='w', transform=ax.transAxes))
            ax.text(.91, .7, 'PD', transform=ax.transAxes, va='center', ha='left')

        
        if i == 5:
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_1_coord),color='k',transform=ax.transAxes))
            ax.add_patch(plt.Polygon(((parent_circle_coord),child_2_coord),color='k',transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(parent_circle_coord,ellipse_width,ellipse_height, color='k', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_1_coord,ellipse_width,ellipse_height, edgecolor='k', facecolor='w', transform=ax.transAxes))
            ax.add_patch(patches.Ellipse(child_2_coord,ellipse_width,ellipse_height, edgecolor='k', facecolor='w', transform=ax.transAxes))
            ax.text(.91, .7, 'DD', transform=ax.transAxes, va='center', ha='left')

        
        if i == 6:
            ax.text(.80,.85,'CMZ', transform=ax.transAxes, fontweight='bold', va='top')

    if legend_objects_2 == None:
        fig_to_ann.legend((legend_objects_1[1][0],legend_objects_1[0][0]),('Observations',legend_objects_1_name), bbox_to_anchor=(.59,.10,.1,.1), loc='upper left')
    else:
        fig_to_ann.legend((legend_objects_1[1][0],legend_objects_1[0][0],legend_objects_2[0][0]),('Observations',legend_objects_1_name, legend_objects_2_name), bbox_to_anchor=(.59,.1,.1,.1), loc='upper left')
    
    #data group box annotations
    fig_to_ann.patches.extend([plt.Rectangle((0.005,0.41),0.99,0.58,
                                  fill=False, color='k', linestyle='-',
                                  transform=fig_to_ann.transFigure, figure=fig_to_ann)])
    fig_to_ann.patches.extend([plt.Polygon(((0.005, .4),(.995,.4), (.995,.21), (.5,.21), (.5, .01), (.005, .01)),
                                  fill=False, color='k', linestyle=':',
                                  transform=fig_to_ann.transFigure, figure=fig_to_ann)])
    
    fig_to_ann.text(.605, .075, 'Training Data', fontsize=12, transform = fig_to_ann.transFigure, figure=fig_to_ann)
    fig_to_ann.patches.extend([plt.Rectangle((0.6,0.07),0.143,0.029,
                                  fill=False, color='k', linestyle='-',
                                  transform=fig_to_ann.transFigure, figure=fig_to_ann)])
    
    
    fig_to_ann.text(.605, .025, 'Test Data', fontsize=12, transform = fig_to_ann.transFigure, figure=fig_to_ann)
    fig_to_ann.patches.extend([plt.Rectangle((0.6,0.02),0.111,0.029,
                                  fill=False, color='k', linestyle=':',
                                  transform=fig_to_ann.transFigure, figure=fig_to_ann)])
    
    fig_to_ann.subplots_adjust(wspace=0.4, hspace=0.3)
    fig_to_ann.tight_layout()

def AIC_table (filename, o_training_AIC, o_test_AIC, s_training_AIC, s_test_AIC, d_training_AIC, d_test_AIC):
    #writes LaTeX tabular snippet with AIC values
    table_file = open(filename, 'w')
    table_file.write(r'\begin{tabular}{|l|l|l|}'+'\n')
    table_file.write(r'\hline'+'\n')
    table_file.write(r'{\bf Model} & {\bf Training AIC} & {\bf Test AIC} \\ \thickhline'+'\n')
    table_file.write(r'Stochastic (He fit) & $'+f'{o_training_AIC:.2f}'+r'$ & $' + f'{o_test_AIC:.2f}' + r'$\\ \hline'+'\n')
    table_file.write(r'Stochastic (SPSA fit) & $'+f'{s_training_AIC:.2f}'+r'$ & $' + f'{s_test_AIC:.2f}' + r'$\\ \hline'+'\n')
    table_file.write(r'Deterministic (SPSA fit) & $'+f'{d_training_AIC:.2f}'+r'$ & $' + f'{d_test_AIC:.2f}' + r'$\\ \hline'+'\n')
    table_file.write(r'\end{tabular}'+'\n')

if __name__ == "__main__":
    main()
