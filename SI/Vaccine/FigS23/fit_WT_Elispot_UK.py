'''
Fitting the immunogenicity distribution and the SARS-CoV-2 protection efficacy of CoronaVac and BNT162b2 (1 or 2 doses)

Author: Dianjie Li
Built Date: 2022/05/28

Edited:
    (1) 2023/3/30 Add boostrapping for quantifying prediction uncertainty
    (2) 2023/05/22 Compute 95% CI with percentiles method
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

from _functions import *
from _boostrap_functions import *

set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']

# -------------- Load Data -----------------------------
filename = './Elispot_Pfizer_UK.xlsx'
vax=['Pfizer','Pfizer-1dose','CoronaVac']

eff_phase3 = []
eff_phase3_CI = []
for i in range(len(vax)):
    df_tmp = pd.read_excel(filename, sheet_name=vax[i])
    eff_phase3.append(df_tmp['Efficacy'][0])
    eff_phase3_CI.append(df_tmp['Efficacy'][1:3].tolist())
eff_phase3 = np.array(eff_phase3)
eff_phase3_CI = np.array(eff_phase3_CI)


# -------------- Preprocessing Data -------------------------------------
"""
norm T -> normalized specific T level = (IFNg - IFNg0_mean)/IFNg0_mean
norm Ab -> Neutralising Antibody/(Antibody of convalescent)
"""

dt = {'normT': [], 'normAb': [], 'name': [], 'logAb': [], 'logT': []}
# Calculate reference CD8T cell level
CD8_ref=[]
for i in range(len(vax)):
    df_tmp = pd.read_excel(filename, sheet_name=vax[i])
    CD8_ref_list=df_tmp['CD8T_ref']
    CD8_ref_list=CD8_ref_list.loc[~(CD8_ref_list==0)]
    CD8_ref_tmp=10 ** np.nanmean(np.log10(np.array(CD8_ref_list)))
    CD8_ref.append(CD8_ref_tmp)
substract=[1,1,1]
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
limit = [[20, 0],[20, 0],[0, 0]]
for i in range(len(dt['name'])):
    T_tmp = dt['normT'][i]
    limit[i][1] = np.min(T_tmp[T_tmp > 0])
    if vax[i]=='CoronaVac':
        Ab_tmp = dt['normAb'][i]
        limit[i][0] = np.min(Ab_tmp[Ab_tmp > 0])
    else:
        limit[i][0] = limit[i][0]/GMT_Ab[i]
dt.update({'limit':limit})

# -------- Fit distribution and e threshold -------------------------------
mu, sigma, dt = fit_distribution(dt)

eth_fit = fit_e_threshold(eff_phase3, mu, sigma, e0=-1)

# eth_fit = fit_e_threshold(eff_phase3[0:2], mu[0:2], sigma[0:2], e0=-1)

eff_pred = np.array([efficacy(mu[i], sigma[i], eth_fit) for i in range(len(mu))])

# -------- Bootstrapping for prediction uncertainty -----------------------
btstrp = boostrap_method(dt, eth_fit)
eff_boostrap = btstrp.boostrap(N_boostrap = 500)

eff_std = np.array([np.std(eff_boostrap[:,k]) for k in range(len(eff_phase3))])
eff_mean = np.array([np.mean(eff_boostrap[:,k]) for k in range(len(eff_phase3))])
eff_CI = np.array([np.percentile(eff_boostrap[:,k], [5, 95]) for k in range(len(eff_phase3))])

plt.figure()
for k in range(3):
    plt.hist(eff_boostrap[:,k], 
             label=dt['name'][k]+' (%.2f, %.2f)' %(eff_CI[k,0],eff_CI[k,1]))
plt.legend()

# -------- Save prediction results ----------------------------------------
pred_save = {'Vaccine':dt['name'],
             'Efficacy':eff_phase3,
             'CI_up':eff_phase3_CI[:,1],
             'CI_down':eff_phase3_CI[:,0],
             'Prediction':eff_pred,
             'Mean_pred':eff_mean,
             'Std_pred':eff_std,
             'CI_up_pred':eff_CI[:,1],
             'CI_down_pred':eff_CI[:,0]
             }
pred_save = pd.DataFrame.from_dict(pred_save)
# pred_save.to_excel('Pred_WT_Elispot.xlsx')

# -------- Draw fitting results -------------------------------------------
dt_joint = {'x': [], 'y': [], 'Vaccine (Efficacy)': []}
splsize = 10000

# Prediction vs Vaccine trials
plt.figure(figsize=(4,2.5))
plt.plot([0,100],[0,100],'--',color='gray')
for i in range(len(eff_pred)):
    yerr_tmp = np.abs(np.array([[eff_phase3_CI[i,0],],[eff_phase3_CI[i,1],]])-eff_phase3[i])
    xerr_tmp = np.array([[eff_pred[i]-eff_CI[i,0]],[-eff_pred[i]+eff_CI[i,1]]])
    plt.errorbar(x=eff_pred[i]*100,y=eff_phase3[i]*100,
                 yerr=yerr_tmp*100,xerr=xerr_tmp*100,
                 marker = 'o',
                 label = vax[i])
plt.legend(loc='lower left',bbox_to_anchor=(1,0.1))
plt.xlabel('Prediction (%)')
plt.ylabel('Efficacy from trials (%)')
plt.xlim([0,105])
plt.ylim([0,105])
plt.grid()
plt.tight_layout()

# sampling data of fitted distribution
for i in range(len(dt['name'])):
    xtmp, ytmp = two_norm_sample(mu[i], sigma[i], splsize)
    dt_joint['x'] = dt_joint['x']+xtmp.tolist()
    dt_joint['y'] = dt_joint['y']+ytmp.tolist()
    dt_joint['Vaccine (Efficacy)'] = \
        dt_joint['Vaccine (Efficacy)'] \
        + [dt['name'][i]+' (%.2f%%)' % (eff_pred[i]*100)]*splsize
    if i == 0:
        eff1 = eff_pred[i] * 100
    if i == 1:
        eff2 = eff_pred[i] * 100

# protection line
g = sns.jointplot(data=dt_joint, x='x', y='y',
                  hue='Vaccine (Efficacy)', kind="kde", height=3.6, space=0,
                  alpha=0.5,
                  levels=6, fill=True, palette=sns.set_palette(set2colors))
ax = g.ax_joint
xtmp = np.linspace(-5, 5, 100)
ytmp = eff_pred-xtmp
ax.plot(xtmp, ytmp, 'k--', lw=2.5, alpha=0.8, zorder=100, label=eff_pred)
ax.text(x=0.5, y=eff_pred-0.5+0.2,
        s='Protection'+'\n'+'Border',
        fontsize=14,
        rotation=-25, rotation_mode="anchor")
print(dt['name'])

# Plot original individual-level data
for i in range(len(mu)):
    xtmp = np.copy(dt['logAb'][i])
    ytmp = np.copy(dt['logT'][i])
    if len(xtmp) != len(ytmp):
        len_sample = min([len(xtmp), len(ytmp)])
        np.random.shuffle(xtmp)
        np.random.shuffle(ytmp)
        xtmp = np.copy(xtmp[:len_sample])
        ytmp = np.copy(ytmp[:len_sample])
    else:
        np.random.shuffle(xtmp)

    for xx, yy in zip(xtmp, ytmp):
        if xx+yy >= eff_pred:
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
# fig.savefig('efficacy_WT.png')
# fig.savefig('efficacy_WT.svg')

plt.show()
