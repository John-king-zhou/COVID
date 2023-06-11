'''
Fitting dynamic data of the immunogenicity distribution and the SARS-CoV-2 protection efficacy against 3 variants.

Author: Dianjie Li
Built Date: 2022/05/30

Edited:
    (1) 2023/3/31 Add boostrapping for quantifying prediction uncertainty
    (2) 2023/05/22 Compute 95% CI with percentiles method
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from _functions import *
from _boostrap_functions import *

set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']

# -------------- Load Data -----------------------------
filename = './Variants_Elispot.xlsx'
vax=['Pfizer','ChAdOx1']

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
substract = [0,1]
CD8_ref[0] = CD8_ref[1]

#load data
GMT_Ab=np.zeros(len(vax))
for i in range(len(vax)):
    df = pd.read_excel(filename, sheet_name=vax[i])
    dt['normT'].append(np.array(df['CD8T']/CD8_ref[i]-substract[i]))
    GMT_Ab[i] = 10 ** np.nanmean(np.log10(np.array(df['ConvAb'])))
    dt['normAb'].append(np.array(df['Ab']) / GMT_Ab[i])
    dt['name'].append(vax[i])

# Limitations of detection for antibody titers and T cell levels
limit = [[20, 0],[1, 0]]
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

# fit distribution of Ab and T
mu, sigma, dt = fit_distribution(dt)

# load vaccine efficacy
variant_seq = ['Alpha', 'Delta', 'Omicron']
set_mrk = ['^', 's', 'o']
eth_fit = []
eff_phase3 = []
eff_phase3_CI = []

for j in range(len(variant_seq)):
    variant = variant_seq[j]
    eff_trials = []
    eff_trials_CI = []
    for i in range(len(vax)):
        df_tmp = pd.read_excel(filename, sheet_name=vax[i])
        eff_trials.append(df_tmp[variant][0])
        eff_trials_CI.append(df_tmp[variant][1:3].tolist())
    eff_phase3.append(eff_trials)
    eff_phase3_CI.append(eff_trials_CI)

    # fit efficacy
    eth_tmp = fit_e_threshold(eff_phase3[-1], mu, sigma, e0=1)
    eth_fit.append(eth_tmp)


eff_std = []
eff_pred = []
eff_mean = []
eff_CI = []
for j in range(len(variant_seq)):
    print('\n'+variant_seq[j])
    eff_pred.append([efficacy(mu[i], sigma[i], eth_fit[j]) for i in range(len(mu))])
    # boostrap
    btstrp = boostrap_method(dt, eth_fit[j])
    eff_boostrap = btstrp.boostrap(N_boostrap = 500)
    eff_std.append(np.array([np.std(eff_boostrap[:,k]) for k in range(len(vax))]))
    eff_mean.append(np.array([np.mean(eff_boostrap[:,k]) for k in range(len(vax))]))
    eff_CI.append(np.array([np.percentile(eff_boostrap[:,k], [5, 95]) for k in range(len(vax))]))

    plt.figure()
    for k in range(2):
        plt.hist(eff_boostrap[:,k], 
                label=dt['name'][k]+' (%.2f, %.2f)' %(eff_CI[-1][k,0],eff_CI[-1][k,1]))
    plt.legend()


# -------- Save prediction results ----------------------------------------
pred_save = {'Vaccine':[],'Strain':[],
             'Efficacy':[],'CI_up':[],'CI_down':[],
             'Prediction':[],'Mean_pred':[],'Std_pred':[],
             'CI_up_pred':[],'CI_down_pred':[],
             }
for i,var_tmp in enumerate(variant_seq):
    for j,vax_tmp in enumerate(dt['name']):
        pred_save['Vaccine'].append(vax_tmp)
        pred_save['Strain'].append(var_tmp)
        pred_save['Efficacy'].append(eff_phase3[i][j])
        pred_save['CI_up'].append(eff_phase3_CI[i][j][1])
        pred_save['CI_down'].append(eff_phase3_CI[i][j][0])
        pred_save['Prediction'].append(eff_pred[i][j])
        pred_save['Mean_pred'].append(eff_mean[i][j])
        pred_save['Std_pred'].append(eff_std[i][j])
        pred_save['CI_up_pred'].append(eff_CI[i][j][1])
        pred_save['CI_down_pred'].append(eff_CI[i][j][0])

pred_save = pd.DataFrame.from_dict(pred_save)
# pred_save.to_excel('Variants_prediction.xlsx')

# ---------------   Figures  ----------------------------------

# sampling data of fitted distribution
dt_joint = {'x': [], 'y': [], 'Vaccine (Efficacy)': []}
splsize = 10000
for i in range(len(dt['name'])):
    xtmp, ytmp = two_norm_sample(mu[i], sigma[i], splsize)
    dt_joint['x'] = dt_joint['x']+xtmp.tolist()
    dt_joint['y'] = dt_joint['y']+ytmp.tolist()
    dt_joint['Vaccine (Efficacy)'] = \
        dt_joint['Vaccine (Efficacy)'] \
        + [dt['name'][i]]*splsize

# protection line
g = sns.jointplot(data=dt_joint, x='x', y='y',
                  hue='Vaccine (Efficacy)', kind="kde", height=3.6, space=0,
                  alpha=0.5,
                  levels=6, fill=True, palette=sns.set_palette(set2colors))
ax = g.ax_joint
for i in range(len(eth_fit)):
    eth_tmp=eth_fit[i]
    xtmp = np.linspace(-5, 5, 20)
    ytmp = eth_tmp-xtmp
    ax.plot(xtmp, ytmp, ls='--', color='k',
            marker=set_mrk[i], markersize=5,
            lw=1.5, alpha=0.8, zorder=100, 
            label=variant_seq[i])
# ax.text(x=0.5, y=eth_tmp-0.5+0.2,
#         s='Protection'+'\n'+'Border',
#         fontsize=14,
#         rotation=-25, rotation_mode="anchor")
print(dt['name'])

# Plot original individual-level data
# for i in range(len(mu)):
#     xtmp = np.copy(dt['logAb'][i])
#     ytmp = np.copy(dt['logT'][i])
#     if len(xtmp) != len(ytmp):
#         len_sample = min([len(xtmp), len(ytmp)])
#         np.random.shuffle(xtmp)
#         np.random.shuffle(ytmp)
#         xtmp = np.copy(xtmp[:len_sample])
#         ytmp = np.copy(ytmp[:len_sample])
#     else:
#         np.random.shuffle(xtmp)

#     for xx, yy in zip(xtmp, ytmp):
#         if xx+yy >= eth_tmp:
#             ax.scatter(xx, yy, s=40, marker='o',
#                        facecolor=set2colors[i], edgecolor='w')
#         else:
#             ax.scatter(xx, yy, s=80, marker='X', facecolor='w')
#             ax.scatter(xx, yy, s=40, marker='x', facecolor=set2colors[i])


ax.set_xlim([-2, 2])
ax.set_xticks([-2, -1, 0, 1, 2])
ax.set_ylim([-1.5, 2])
ax.set_yticks([ -1, 0, 1, 2])
ax.set_xlabel('$log_{10}$Ab (Normalized)')
ax.set_ylabel('$log_{10}$T (Normalized)')
fig = g.fig
fig.set_size_inches(4, 4)

# to save figures
# fig.savefig('efficacy_variants.png')
# fig.savefig('efficacy_variants.svg')

plt.show()
