import os
import numpy as np
import matplotlib.pyplot as plt
from fileinput import filename
from scipy.stats import lognorm
from scipy.stats import gamma
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

#Model parameters
gomes_normal_sigma = .32839
gomes_normal_mu = 3.9716

he_gamma_shape = 2
he_gamma_scale = 1
he_gamma_offset = 4

def main():
    gomes_x_range = np.arange(1,151,1)
    gomes_y = lognorm.pdf(gomes_x_range, gomes_normal_sigma, 0, np.exp(gomes_normal_mu))
    
    he_x_range = np.arange(1,16,.1)
    he_y = gamma.pdf(he_x_range, he_gamma_shape, he_gamma_offset, he_gamma_scale)
    
    fig, (ax_g, ax_h) = plt.subplots(1,2,figsize=(6,2))
    
    ax_g.plot(gomes_x_range,gomes_y,'k-')
    ax_h.plot(he_x_range, he_y, 'k-')
    
    plt.sca(ax_g)    
    plt.xlabel("Cell cycle duration (h)")
    plt.ylabel("Probability")
    
    plt.sca(ax_h)
    plt.xlabel("Cell cycle duration (h)")
    plt.ylabel("Probability")
    
    plt.tight_layout()
    plt.savefig('/home/main/Desktop/utility.png', transparent=True)
    plt.show()
    
if __name__ == "__main__":
    main()
