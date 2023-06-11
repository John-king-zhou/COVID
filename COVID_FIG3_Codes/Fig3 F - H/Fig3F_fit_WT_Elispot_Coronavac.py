'''
Fitting the immunogenicity distribution and the SARS-CoV-2 protection efficacy of Sputinik V, CoronaVac and BNT162b2

Author: Zhengqing Zhou
Built Date: 2022/05/28
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm
from scipy.stats.mstats import gmean
from scipy.optimize import minimize
import seaborn as sns
from scipy.integrate import quad
from scipy.optimize import curve_fit

from _functions import *

set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']

# -------------- Load Data -----------------------------
filename = './Elispot.xlsx'
vax=['CoronaVac','Pfizer']
eff_phase3 = np.array([0.5065, 0.95,]) # Efficacies from clinical trials


# -------------- Preprocessing Data -------------------------------------
"""
norm T -> normalized specific T level = (IFNg - IFNg0_mean)/IFNg0_mean
norm Ab -> Neutralising Antibody/(Antibody of convalescent)
"""

dt = {'normT': [], 'normAb': [], 'name': [], 'logAb': [], 'logT': []}
# Calculate reference CD8T cell level
df_CoronaVac = pd.read_excel(filename, sheet_name=vax[0])
CD8_ref_list=df_CoronaVac['CD8T_ref']
CD8_ref_list=CD8_ref_list.loc[~(CD8_ref_list==0)]
CD8_ref_CoronaVac=10 ** np.nanmean(np.log10(np.array(CD8_ref_list)))
df_ChAdOx1 = pd.read_excel(filename, sheet_name='ChAdOx1')
CD8_ref_list=df_ChAdOx1['CD8T_ref']
CD8_ref_list=CD8_ref_list.loc[~(CD8_ref_list==0)]
CD8_ref_ChAdOx1=10 ** np.nanmean(np.log10(np.array(CD8_ref_list)))

CD8_ref=[CD8_ref_CoronaVac,CD8_ref_ChAdOx1]
substract=[1,0]
print('CD8',CD8_ref)
#load data
GMT_Ab=np.zeros(len(vax))
for i in range(len(vax)):
    df = pd.read_excel(filename, sheet_name=vax[i])
    dt['normT'].append(np.array(df['CD8T']/CD8_ref[i]-substract[i]))
    GMT_Ab[i] = 10 ** np.nanmean(np.log10(np.array(df['ConvAb'])))
    dt['normAb'].append(np.array(df['Ab']) / GMT_Ab[i])
    dt['name'].append(vax[i])

# Limitations of detection for antibody titers and T cell levels
limit = [[0, 0], [20, 0]]
for i in range(len(dt['name'])):
    T_tmp = dt['normT'][i]
    limit[i][1] = np.min(T_tmp[T_tmp > 0])
    if limit[i][0]==0:
        Ab_tmp = dt['normAb'][i]
        limit[i][0] = np.min(Ab_tmp[Ab_tmp > 0])
    else:
        limit[i][0] = limit[i][0]/GMT_Ab[i]
limit[1][1] = limit[0][1]  # assuming same limit for coronavac and pfizer
dt.update({'limit':limit})


# -------- Fit distribution and e threshold -------------------------------
mu, sigma, dt = fit_distribution(dt)

eth_fit = fit_e_threshold(eff_phase3, mu, sigma)

# -------- Draw fitting results -------------------------------------------
dt_joint = {'x': [], 'y': [], 'Vaccine (Efficacy)': []}
splsize = 10000
eth_tmp = eth_fit
eff_tmp = [efficacy(mu[i], sigma[i], eth_tmp) for i in range(len(mu))]

# sampling data of fitted distribution
for i in range(len(dt['name'])):
    xtmp, ytmp = two_norm_sample(mu[i], sigma[i], splsize)
    dt_joint['x'] = dt_joint['x']+xtmp.tolist()
    dt_joint['y'] = dt_joint['y']+ytmp.tolist()
    dt_joint['Vaccine (Efficacy)'] = \
        dt_joint['Vaccine (Efficacy)'] \
        + [dt['name'][i]+' (%.2f%%)' % (eff_tmp[i]*100)]*splsize
    if i == 0:
        eff1 = eff_tmp[i] * 100
    if i == 1:
        eff2 = eff_tmp[i] * 100

# protection line
g = sns.jointplot(data=dt_joint, x='x', y='y',
                  hue='Vaccine (Efficacy)', kind="kde", height=3.6, space=0,
                  levels=6, fill=True, palette=sns.set_palette(set2colors))
ax = g.ax_joint
xtmp = np.linspace(-5, 5, 100)
ytmp = eth_tmp-xtmp
ax.plot(xtmp, ytmp, 'k--', lw=2.5, alpha=0.8, zorder=100, label=eff_tmp)
ax.text(x=0.5, y=eth_tmp-0.5+0.2,
        s='Protection'+'\n'+'Border',
        fontsize=14,
        rotation=-25, rotation_mode="anchor")
print(dt['name'])
# Plot original individual-level data
for i in range(len(mu)):
    xtmp = np.copy(dt['logAb'][i])
    ytmp = np.copy(dt['logT'][i])
    if dt['name'][i] != 'CoronaVac':
        np.random.shuffle(xtmp)
    elif len(xtmp) != len(ytmp):
        xtmp = xtmp[:len(ytmp)]

    for xx, yy in zip(xtmp, ytmp):
        if xx+yy >= eth_tmp:
            ax.scatter(xx, yy, s=40, marker='o',
                       facecolor=set2colors[i], edgecolor='w')
        else:
            ax.scatter(xx, yy, s=80, marker='X', facecolor='w')
            ax.scatter(xx, yy, s=40, marker='x', facecolor=set2colors[i])
ax.set_xlim([-2.5, 1.8])
ax.set_xticks([-2, -1, 0, 1])
ax.set_ylim([-3.2, 5.8])
ax.set_xlabel('$log_{10}$Ab (Normalized)')
ax.set_ylabel('$log_{10}$T (Normalized)')
fig = g.fig
fig.set_size_inches(4, 4)

# to save figures
# fig.savefig('efficacy.png')
# fig.savefig('efficacy.svg')

plt.show()
