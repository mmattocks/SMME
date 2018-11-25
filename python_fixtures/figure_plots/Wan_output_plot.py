import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from PIL import Image
from io import BytesIO

#PLoS formatting stuff
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

def main():
    wan_output_dir = '/home/main/git/chaste/projects/ISP/python_fixtures/testoutput/WanOutput'
    wan_results_list = []
    for root,dirs,files in os.walk(wan_output_dir):
        for name in files:
            if name == 'celltypes.dat':
                results = np.loadtxt(os.path.join(root,name), usecols=2)
                if len(results) == 8569:
                    wan_results_list.append(results)
      
    wan_sim_results = np.array(wan_results_list)
    wan_sim_mean = np.mean(wan_sim_results,axis=0)
    wan_sim_95CI = 2*np.std(wan_sim_results,axis=0)
    wan_sim_x_seq = np.arange(72,8641,1)/24
     
    empirical_pop_mean = np.array([792.0, 768.4, 906.1, 1159.7, 1630.0, 3157.9, 3480.2, 4105.1, 1003.0, 477.2, 438.8088611111])
    empirical_pop_95CI = 2*np.array([160.1, 200.1, 244.5, 477.6, 444.3, 1414.3, 472.1, 1169.7, 422.8, 367.5, 294.8])
    empirical_x_seq = np.array([72,120,192,288,408,552,720,1440,2160,4320,8640])/24
    
    fig, ax = plt.subplots(1,1)
    
    empirical = plt.plot(empirical_x_seq, empirical_pop_mean, 'g-')
    plt.fill_between(empirical_x_seq, (empirical_pop_mean - empirical_pop_95CI), (empirical_pop_mean + empirical_pop_95CI), alpha=0.1, edgecolor='#008000', facecolor='#00FF00')
     
    wan_simulator = plt.plot(wan_sim_x_seq, wan_sim_mean, 'm-')
    plt.fill_between(wan_sim_x_seq, (wan_sim_mean - wan_sim_95CI), (wan_sim_mean + wan_sim_95CI), alpha=0.1, edgecolor='#800080', facecolor='#FF00FF')
    
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    
    plt.ylabel("CMZ population size")
    plt.xlabel("Retina age (dpf)")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(40))
    ax.legend((empirical[0], wan_simulator[0]), ('Observed CMZ population', '"Wan-type" simulated CMZ population'))
    
    png_memory = BytesIO()
    fig.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('wan_fig.tiff')
    png_memory.close()
        
    plt.show()

if __name__ == "__main__":
    main()