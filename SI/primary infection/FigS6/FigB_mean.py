'''

Compare Mean VL and Immune data between clinical data and simulations.

Lineplot and errorbarplot

Date: 2023/05/01

'''
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from Read_data import Data

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
mode_list = [2,3]


# load clinical data
data = Data()

dt_mean = {'group':[], 'Imm_val':[], 'Imm_idx':[],
           'weeks':[], 'days':[]}

for id in data.dataset['id']:
    dt_tmp = data.single_patient_data(id)
    if 'Mild' in dt_tmp['group']:
        grp = 'Mild'
    else:
        grp = 'Moderate/Severe'

    for t_idx, idx in zip(['t_VL','t_Imm','t_Imm', 't_Imm'], ['lgVL','Ab','T cell', 'Ab x T']):
        for t, val in zip(dt_tmp[t_idx], dt_tmp[idx]):
            dt_mean['weeks'].append(int(t//7)+1)
            dt_mean['days'].append(t)
            dt_mean['Imm_idx'].append(idx)
            dt_mean['Imm_val'].append(val)
            dt_mean['group'].append(grp)


# plot clinical data
fig_ax = []
clr_list = [ggcolors[i] for i in mode_list]
ylim_list = [[1, 6], [-2, 102], [-100,3000], [-100,220000]]
y_label = [r'$log_{10}$(VL)', '% inhibition', 
           r'Total SFU/$10^6$ PBMC', 'Ab X T']
title_seq = ['Viral load', 'Humoral response', 
             'Cellular response', 'Ab X T']

fig, ax_all = plt.subplots(4,2,figsize=(6,12))
for i_idx, idx in enumerate(['lgVL','Ab','T cell','Ab x T']):
    df = pd.DataFrame.from_dict(dt_mean)
    df = df[df['Imm_idx']==idx]
    if idx == 'lgVL':
        df = df[df['weeks']<5]

    ax = ax_all[i_idx,:]
    for i, grp in zip(range(2), ['Mild', 'Moderate/Severe']):
        df_tmp = df[df['group']==grp]
        sns.pointplot(data=df_tmp, x="weeks", y='Imm_val', 
                      capsize=.2, color=clr_list[i], 
                      errorbar="sd", errwidth=1,
                      scale=0.7, legend=None,
                      linestyles='', ax=ax[i],)
        ax[i].errorbar(x=1000,y=1000,yerr=10,capsize=3,color=clr_list[i],
                       label=grp+' patients',marker='o',ls='')

        ax[i].set_ylim(ylim_list[i_idx])
        ax[i].set_ylabel(y_label[i_idx])
        ax[i].set_title(title_seq[i_idx])
        ax[i].set_xlabel('Days post symptom onset')
        ax[i].set_xlim([-1.2, 3.2])
        ax[i].set_xticks(np.arange(-1,4,1)-3/7, np.arange(-1,4,1)*7)

    fig_ax.append([fig, ax])



# Load simulatino results

# model --> clinical 
def model_to_clinical(mean_x, std_x, type):
    if type == 0:  # lgVL
        lgVL_add = 2
        new_lgVL = mean_x[0]+lgVL_add
        # new_lgVL[new_lgVL<1] = 1
        return new_lgVL, std_x[0]
    elif type == 1:  # Ab
        Ab_scale = 10
        return mean_x[30]/Ab_scale, std_x[30]/Ab_scale
    elif type == 2:  # IFNg
        IFNg_scale = 1500
        return (mean_x[19] + mean_x[20])*IFNg_scale, \
                (std_x[19] + std_x[20])*IFNg_scale
    

time=np.arange(0, 50, 0.1)
time_symptom_onset = (time-10)/7
for i in range(3):
    fig, ax = fig_ax[i]
    for j, mode in enumerate(mode_list):

        Aver_Results=np.loadtxt('Mean%i.txt'%(mode))
        Std_Results=np.loadtxt('Std%i.txt'%(mode))

        Imm, Imm_Std = model_to_clinical(Aver_Results, Std_Results, i)
        
        ax[j].plot(time_symptom_onset, Imm, c=ggcolors[mode],
                   alpha=0.5, label='Mode%i (Averaged)' %mode)
        ax[j].fill_between(time_symptom_onset, Imm-Imm_Std, Imm+Imm_Std, 
                           facecolor=ggcolors[mode], alpha=0.1)
        if i==0:
            ax[j].legend()

fig.subplots_adjust(top=0.912,
                    bottom=0.22,
                    left=0.108,
                    right=0.979,
                    hspace=0.6,
                    wspace=0.35)

# fig.savefig('Fig_Compare_2.svg', dpi=200)
# fig.savefig('Fig_Compare_2.png', dpi=200)
        
# e vs Ab x T

