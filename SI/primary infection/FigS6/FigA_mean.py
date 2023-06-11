'''

Mean VL and Immune data of Tan 2021.

Date: 2023/04/30
Edited: 2023/05/13 three subfigure together with legends

'''
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from Read_data import Data

data = Data()

dt_mean = {'group':[], 'Imm_val':[], 'Imm_idx':[],
           'weeks':[]}

for id in data.dataset['id']:
    dt_tmp = data.single_patient_data(id)
    if 'Mild' in dt_tmp['group']:
        grp = 'Mild'
    else:
        grp = 'Moderate/Severe'

    for t_idx, idx in zip(['t_VL','t_Imm','t_Imm'], ['lgVL','Ab','T cell']):
        for t, val in zip(dt_tmp[t_idx], dt_tmp[idx]):
            dt_mean['weeks'].append(int(t//7)+1)
            dt_mean['Imm_idx'].append(idx)
            dt_mean['Imm_val'].append(val)
            dt_mean['group'].append(grp)


# Figure

fig, ax = plt.subplots(3,1,figsize=(4,7))
y_label = [r'$log_{10}$(VL)', '% inhibition', 
           r'Total SFU/$10^6$ PBMC']

for i_idx, idx in enumerate(['lgVL','Ab','T cell']):
    df = pd.DataFrame.from_dict(dt_mean)
    df = df[df['Imm_idx']==idx]
    df = df[df['weeks']<10]

    sns.barplot(data=df, x='weeks', y='Imm_val', hue='group',
                capsize=.2, lw=2, edgecolor=".5",
                palette='RdBu_r',ax=ax[i_idx])
    sns.stripplot(x="weeks", y="Imm_val", hue='group', data=df,
                  palette='RdBu_r', size=5, linewidth=1,
                  dodge=True, alpha=.8, zorder=100, ax=ax[i_idx])
    # ax.legend([],[],frameon=False)

    ax[i_idx].set_ylim([0, max(df['Imm_val'])*1.1])
    ax[i_idx].set_ylabel(y_label[i_idx])
    ax[i_idx].set_xlabel('Weeks post symptom onset')
    ax[i_idx].set_xlim([-0.5, 3.5])
fig.tight_layout()

# fig.savefig('Fig_Compare_1.svg', dpi=200)
# fig.savefig('Fig_Compare_1.png', dpi=200)

