import os
import numpy as np
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from PIL import Image
from io import BytesIO


#PLoS formatting stuff
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

def main():
    #Read & parse raw data file
    raw_data = pd.read_excel('/home/main/git/chaste/projects/ISP/empirical_data/cumulative_edu.xlsx')
    #group data by dorsal/ventral
    dv_grouped = raw_data.groupby('D/V')
    #make new dataframe for totals and labelled fraction
    totals_frame = dv_grouped.get_group('D')
    #reset indices to ventral frame for addition operatinos
    totals_frame = totals_frame.set_index(dv_grouped.get_group('V').index)
    #sum D + V PCNA & EdU counts
    totals_frame.PCNA = totals_frame.PCNA + dv_grouped.get_group('V').PCNA
    totals_frame.EdU = totals_frame.EdU + dv_grouped.get_group('V').EdU
    #create new column for labelled fraction
    totals_frame['labelled_fraction'] = totals_frame.EdU / totals_frame.PCNA
    
    #setup x for linear regression with y-intercept constant
    X = totals_frame.Time
    X = sm.add_constant(X)
    
    model = sm.OLS(totals_frame.labelled_fraction, X).fit()
    print(model.summary())
    fit_results = model.params
    
    tc = 1/fit_results[1]
    ts = tc * fit_results[0]
    
    sns.regplot(totals_frame.Time,totals_frame.labelled_fraction)

    plt.text(7,.2,"Tc: " + f'{tc:.2f}')
    plt.text(7,.1,"Ts: " + f'{ts:.2f}')

    plt.xlabel("Time (h)")
    plt.ylabel("Fraction of CMZ RPCs labelled by EdU")

    plt.savefig("cumulative_edu.png")
    png_memory = BytesIO()
    plt.savefig(png_memory, format='png', dpi=600)
    PILpng = Image.open(png_memory)
    PILpng.save('cumulative_edu.tiff')
    png_memory.close()
    
    plt.show()

if __name__ == "__main__":
    main()
