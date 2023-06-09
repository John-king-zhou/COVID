#data from Goyal et al. Science adv. 2020
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from FIT1 import FIT2, Imm1, Imm2
import warnings

warnings.filterwarnings('ignore')

set2colors = ['#66c2a5','#fc8d62','#a6d854','#e78ac3','#ffd92f','#e5c494','#b3b3b3']
df=pd.read_csv('Goyal2020.csv')
'''
Individuals
'''
indivs=list(set(df.loc[:,'ID']))
'''
preclude the ones with two or less data points
'''
new_indivs=[]
dfs=[]
for id in range(len(indivs)):
    indiv=indivs[id]
    inds=df.index[df.loc[:,'ID']==indiv].tolist()
    df_i=df.iloc[inds,:]
    time=np.array(df_i.loc[:,'time'])
    if len(time)<3:
        continue
    else:
        new_indivs.append(indiv)
        dfs.append(df_i)
print(len(dfs))
dfs=pd.concat(dfs)
dfs=dfs.reset_index(drop=True)
df2={'ID':[],'e':[],'e_type':[],'lgvmax':[],'recovery':[]}
fig,axes=plt.subplots(5,5,figsize=(6,6))
axes=axes.flat
for i in range(len(new_indivs)):
    ratio = i / len(new_indivs)
    rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
    rat_str = ''.join(rat_str)
    print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
    id=new_indivs[i]
    ax=axes[i]
    inds=dfs.index[dfs.loc[:,'ID']==id].tolist()
    df_i=dfs.iloc[inds,:]
    time=np.array(df_i.loc[:,'time'])
    t0=time[0]
    y = np.array(df_i.loc[:, 'logv'])
    results=FIT2(time,y)
    e1,e2=np.nan,np.nan
    t_rec=0
    if results[0]:#when there is a rising phase
        popt=results[1]
        pcov=results[2]
        Imm3=results[3]
        time=time-time[0]
        e1=popt[1]*1
        e2=popt[2]*3
        time33 = np.arange(0, 50, 0.1)
        y33 = Imm3(time33, *popt)
        if len(time33[y33 == 2]) != 0:
            t_rec = time33[y33 == 2][0]-time33[np.argmax(y33)]
        else:
            t_rec = 50
        ax.plot(time33 + t0, y33, color=set2colors[1])
        ax.scatter(time+t0, y, s=30, marker='x', color=set2colors[1])
        ax.text(x=18,y=5.5,s='$\epsilon$1=%1.1f \n$\epsilon$2=%1.1f' % (e1, e2),fontsize=7,ha='center',va='center')
        ax.text(x=18, y=8.5, s='ID=%s' % (id),fontsize=7,ha='center',va='center')
    else:
        popt=results[2]
        pcov=results[3]
        if results[1] == 1:
            t_rec=0
            time = time - time[0]
            e2 = popt[1] * 1
            time11 = np.arange(0, 50, 0.1)
            y11=Imm1(time11, *popt)
            ax.plot(time11+t0, y11, color=set2colors[1])
        else:
            time = time - time[0]
            e2 = popt[1] * 3
            time22 = np.arange(0, 50, 0.1)
            y22=Imm2(time22, *popt)
            if len(time22[y22==2])!=0:
                t_rec=time22[y22==2][0]
            else:
                t_rec=50
            ax.plot(time22+t0, y22, color=set2colors[1])
        ax.scatter(time+t0, y, s=30, marker='x', color=set2colors[1])
        ax.text(x=18, y=5.5, s='$\epsilon$%i=%1.1f' % (results[1],e2),fontsize=7,ha='center',va='center')
        ax.text(x=18, y=8.5, s='ID=%s' % (id),fontsize=7,ha='center',va='center')
    ax.set_xlim([-2, 25])
    ax.set_xticks([0, 10, 20])
    ax.set_ylim([0.5, 10.4])
    ax.set_yticks([2,5,8])
    if i%5==0:
        ax.set_yticklabels([2,5,8])
    else:
        ax.set_yticklabels([])
    if i>19:
        ax.set_xticklabels([0, 10, 20])
    else:
        ax.set_xticklabels([])
    if not np.isnan(e1):
        df2['ID']+=[id,]
        df2['e']+=[e1,]
        df2['e_type']+=[1,]
        df2['lgvmax']+=[np.max(y),]
        df2['recovery']+=[t_rec,]
    if not np.isnan(e2):
        df2['ID']+=[id,]
        df2['e']+=[e2,]
        df2['e_type']+=[2,]
        df2['lgvmax']+=[np.max(y),]
        df2['recovery']+=[t_rec,]
df2=pd.DataFrame(df2)
df2.to_csv('Goyal_e.csv')
fig.subplots_adjust(bottom=0.1,top=0.98,left=0.1,right=0.98)
fig.text(0.5,0.05,s='Days post symptom onset',ha='center',va='center')
fig.text(0.05,0.5,s='Goyal 2020 Science Advances \n log$_{10}$nCoV',ha='center',va='center',rotation=90)
fig.savefig('Fit_All_Goyal.png')
fig.savefig('Fit_All_Goyal.svg')
plt.show()