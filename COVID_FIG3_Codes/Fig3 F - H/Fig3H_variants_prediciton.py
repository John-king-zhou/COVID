import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sqlalchemy import column
set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
marker_seq = ['^', 's', 'o']

data = pd.read_excel('Variants_prediction.xlsx')
vax = data['Vaccine'].drop_duplicates()
variants = data['Strain'].drop_duplicates()

plt.figure(figsize=(4.8,2.5))
plt.plot([0,200],[0,200],
         '--',color='black',
         lw=1)

for i, v in enumerate(vax):
    for j, var in enumerate(variants):
        data_tmp = data[data['Strain']==var]
        data_tmp = data_tmp[data_tmp['Vaccine']==v].iloc[0,1:].to_dict()

        yerr = np.array([[-data_tmp['CI_down']+data_tmp['Efficacy']], 
                        [data_tmp['CI_up']-data_tmp['Efficacy']]])
        xerr = np.array([[data_tmp['Prediction']-data_tmp['CI_down_pred']],
                         [data_tmp['CI_up_pred']-data_tmp['Prediction']]])

        plt.errorbar(x=data_tmp['Prediction']*100,
                     y=data_tmp['Efficacy']*100,
                     yerr = yerr*100, xerr = xerr*100,
                     ls='',
                     markersize = 6,
                     marker = marker_seq[j],
                     markeredgecolor = 'gray',
                     color = set2colors[i],
                     label = v+' ('+var+')',
                     alpha=0.6)
    
plt.legend(loc='lower left',bbox_to_anchor=(1,0.1))
plt.xlabel('Prediction (%)')
plt.ylabel('Efficacy from trials (%)')
plt.xlim([20,100])
plt.ylim([20,100])
plt.title('Variants')
# plt.grid()
plt.tight_layout()
# plt.savefig('Variants_Prediction.svg')
# plt.savefig('Variants_Prediction.png', dpi=200)