import os
import numpy as np
import matplotlib.pyplot as plt

#AIC & Plotting utility params
event_seeds = 1000
error_samples = 5000 #number of samples to draw when estimating plausibility interval for simulations

def main():
    #LOAD & PARSE KOLMOGOROV COMPLEXITY FILES
    gomes_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/gomes_acss', skiprows=1, usecols=1)
    he_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/he_acss', skiprows=1, usecols=1)
    he_refit_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/heRefit_acss', skiprows=1, usecols=1)
    deterministic_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/deterministic_acss', skiprows=1, usecols=1)
    boije_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/boije_acss', skiprows=1, usecols=1)
    eo_kol = np.loadtxt('/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/EO_acss', skiprows=1, usecols=1)
    
    data_sequence = [gomes_kol,he_kol,he_refit_kol,deterministic_kol,boije_kol,eo_kol]
    
    plt.violinplot(data_sequence, showextrema=True)
    
    plt.ylabel ('Estimated Kolmogorov Complexity')
    
    plt.show()

if __name__ == "__main__":
    main()
