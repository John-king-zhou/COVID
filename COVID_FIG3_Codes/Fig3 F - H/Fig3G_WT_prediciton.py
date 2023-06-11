import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sqlalchemy import column
set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
marker_seq = ['o','s','^','*','p']

data = pd.read_excel('WT_prediction.xlsx')
vax = data['Vaccine']

plt.figure(figsize=(2.8,2.5))
plt.plot([0,200],[0,200],'--',color='black',
         lw=1)

for v,i in zip(vax, range(len(vax))):
    data_tmp = data[data['Vaccine']==v].iloc[0,1:].to_dict()
    yerr = np.array([[-data_tmp['CI_down']+data_tmp['Efficacy']], 
                     [data_tmp['CI_up']-data_tmp['Efficacy']]])
    xerr = np.array([[data_tmp['Prediction']-data_tmp['CI_down_pred']],
                     [data_tmp['CI_up_pred']-data_tmp['Prediction']]])

    plt.errorbar(x=data_tmp['Prediction']*100,
                 y=data_tmp['Efficacy']*100,
                 yerr = yerr*100, xerr = xerr*100,
                 ls='',
                 markersize = 8,
                 marker = marker_seq[i],
                 markeredgecolor = 'gray',
                 color = set2colors[i],
                 label = v)
plt.legend(loc='best',ncol=1,
           )
plt.xlabel('Prediction (%)')
plt.ylabel('Efficacy from trials (%)')
plt.xticks([25,50,75,100])
plt.yticks([25,50,75,100])
plt.xlim([20,102])
plt.ylim([20,102])
plt.title('WT')
# plt.grid()
plt.tight_layout()
plt.savefig('WT_prediction.svg')
plt.savefig('WT_prediction.png',dpi=200)
plt.show()