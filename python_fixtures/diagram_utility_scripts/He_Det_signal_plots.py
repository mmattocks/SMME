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

d_pp_y = [1,0,0,0]
d_pd_y = [0,1,0,0]
d_dd_y = [0,0,1,1]

fig, ((spp,spd,sdd))= plt.subplots(1,3,figsize=(6,2),sharey=True)
fig2, ((dpp,dpd,ddd))= plt.subplots(1,3,figsize=(6,2),sharey=True)

for [pp,pd,dd] in [[spp,spd,sdd],[dpp,dpd,ddd]]:
    for p in [pp,pd,dd]:
        p.set_xlim((0,24))
        p.set_xticks([8,15])
        p.set_ylim((-.1,1.1))
        p.set_yticks(np.arange(0,1.2,.2))
        p.set_xlabel("TiL (h)")
        if p==spp or p==spd or p==sdd:
            p.axvline(8,linestyle='--',color='.1', alpha=.2)
            p.axvline(15,linestyle='--',color='.1', alpha=.2)
        else:
            p.set_xticklabels(["<->", "<->"])
            p.tick_params(axis="x",colors="red")
            p.axvline(8,linestyle='-',color='blue', alpha=.15, linewidth=15)
            p.axvline(15,linestyle='-',color='blue', alpha=.15, linewidth=15)
    pp.set_ylabel("Probability")
    if pp==spp:
        pp.step(he_x,he_pp_y, 'k-', where='post', linewidth=2)
    else:
        pp.step(he_x,d_pp_y, 'r-', where='post', linewidth=2)
    if pd==spd:
        pd.step(he_x,he_pd_y, 'k-', where='post', linewidth=2)
    else:
        pd.step(he_x,d_pd_y, 'r-', where='post', linewidth=2)
    if dd==sdd:
        dd.step(he_x,he_dd_y, 'k-', where='post', linewidth=2)
    else:
        dd.step(he_x,d_dd_y, 'r-', where='post', linewidth=2)

fig.tight_layout()
fig2.tight_layout()
fig.savefig('Hemode.png')
fig2.savefig('Detmode.png')
plt.show()
     