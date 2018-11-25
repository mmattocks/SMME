import os
import numpy as np
import matplotlib.pyplot as plt
from fileinput import filename
from scipy.stats import lognorm
from scipy.stats import gamma
from statsmodels.sandbox.distributions.gof_new import a_st70_upp

#Model parameters
he_x = [0,8,15,25]
he_pp_y = [1,.2,.2,.2]
he_pd_y = [0,.4,0,0]
he_dd_y = [0,.4,.8,.8]

boije_x = [0,3,5,9]
boije_pAtoh7 = [0,.32,0,0]
boije_pPtf1a = [0,.3,0,0]
boije_pNg = [0,0,.8,.8]

def main():
    fig, ((pp,pd,dd))= plt.subplots(1,3,figsize=(6,2),sharey=True)
    
    fig2, ((pAtoh7),(pPtf1a),(pNg)) = plt.subplots(3,1,figsize=(2,4),sharex=True)
    
    pp.set_xlim((0,24))
    pd.set_xlim((0,24))
    dd.set_xlim((0,24))

    pp.set_xticks([8,15])
    pd.set_xticks([8,15])
    dd.set_xticks([8,15])

    pp.set_ylim((-.1,1.1))
    pp.set_yticks(np.arange(0,1.2,.2))
    pp.set_ylabel("Probability")
    pp.set_xlabel("TiL (h)")
    pp.axvline(8,linestyle='--',color='.1', alpha=.2)
    pp.axvline(15,linestyle='--',color='.1', alpha=.2)
    pp.step(he_x,he_pp_y, 'k-', where='post', linewidth=2)

    pd.set_ylim((-.1,1.1))
    pd.set_yticks(np.arange(0,1.2,.2))
    pd.set_xlabel("TiL (h)")
    pd.axvline(8,linestyle='--',color='.1', alpha=.2)
    pd.axvline(15,linestyle='--',color='.1', alpha=.2)
    pd.step(he_x,he_pd_y, 'k-', where='post', linewidth=2)
    
    dd.set_ylim((-.1,1.1))
    dd.set_yticks(np.arange(0,1.2,.2))
    dd.set_xlabel("TiL (h)")
    dd.axvline(8,linestyle='--',color='.1', alpha=.2)
    dd.axvline(15,linestyle='--',color='.1', alpha=.2)
    dd.step(he_x,he_dd_y, 'k-', where='post', linewidth=2)
    
    pAtoh7.set_xlim((0,8))
    pNg.set_xlim((0,8))
    pPtf1a.set_xlim((0,8))

    pAtoh7.set_xticks([3,5])
    pNg.set_xticks([3,5])
    pPtf1a.set_xticks([3,5])

    
    pNg.set_xlabel("Generation")

    
    pAtoh7.set_ylim((-.1,1.1))
    pAtoh7.set_yticks(np.arange(0,1.25,.25))
    pAtoh7.set_ylabel("Probability")
    pAtoh7.axvline(3,linestyle='--',color='.1',alpha=.2)
    pAtoh7.axvline(5,linestyle='--',color='.1',alpha=.2)
    pAtoh7.step(boije_x,boije_pAtoh7, 'k-', where='post', linewidth=2)
    
    pPtf1a.set_ylim((-.1,1.1))
    pPtf1a.set_yticks(np.arange(0,1.25,.25))
    pPtf1a.set_ylabel("Probability")
    pPtf1a.axvline(3,linestyle='--',color='.1',alpha=.2)
    pPtf1a.axvline(5,linestyle='--',color='.1',alpha=.2)
    pPtf1a.step(boije_x,boije_pPtf1a, 'k-', where='post', linewidth=2)
    
    pNg.set_ylim((-.1,1.1))
    pNg.set_yticks(np.arange(0,1.25,.25))
    pNg.set_ylabel("Probability")
    pNg.axvline(3,linestyle='--',color='.1',alpha=.2)
    pNg.axvline(5,linestyle='--',color='.1',alpha=.2)
    pNg.step(boije_x,boije_pNg, 'k-', where='post', linewidth=2)


    plt.tight_layout()
    fig.savefig('Hemode.png', transparent=True)
    fig2.savefig('Boijesignals.png', transparent=True)
    plt.show()
    
if __name__ == "__main__":
    main()
